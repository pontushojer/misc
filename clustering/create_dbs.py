"""
Silmulate randomly generated barcodes.
Based on UMI-Tools pipeline:
https://github.com/CGATOxford/UMI-tools_pipelines/blob/master/notebooks/Simulating_umi_deduping.ipynb
"""

# TODO multiprocessing not working.

import argparse
import random
import collections
import numpy as np
import logging
import time
from tqdm import tqdm
import math
import itertools
import multiprocessing
import dnaio

logger = logging.getLogger(__name__)

counter = collections.Counter()
true_barcodes = set()

WIDTH = 50
MULTIPROCESS_CREATE_MP_LIMIT = 10000

errors_dict = {"A": ["C", "G", "T"],
               "C": ["A", "G", "T"],
               "G": ["C", "A", "T"],
               "T": ["C", "G", "A"]}

translate_dict = {"A": ("A"),
                  "T": ("T"),
                  "C": ("C"),
                  "G": ("G"),
                  "R": ("A", "G"),
                  "Y": ("C", "T"),
                  "M": ("A", "C"),
                  "K": ("G", "T"),
                  "S": ("C", "G"),
                  "H": ("A", "C", "T"),
                  "B": ("C", "G", "T"),
                  "D": ("A", "G", "T"),
                  "V": ("A", "C", "G"),
                  "N": ("A", "G", "C", "T")}


def translate(base):
    return translate_dict[base.upper()]


def get_chunks(iterable, size=10):
    iterator = iter(iterable)
    for first in iterator:
        yield itertools.chain([first], itertools.islice(iterator, size - 1))


def get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--sequence", type=str, default="BDVHBDVHBDVHBDVHBDVH",
                        help="Barcode sequence with degenerate bases. Default: 'BDVHBDVHBDVHBDVHBDVH'")
    parser.add_argument('-n', '--number', type=int, default=100,
                        help="Number of unique barcodes to generate.")
    parser.add_argument('-r', '--reads', type=float, default=None,
                        help="Number of reads to generate. Limits PCR.")
    parser.add_argument('-ep', '--error-rate-pcr', type=float, default=1e-5,
                        help="Error rate in PCR.")
    parser.add_argument('-es', '--error-rate-seq', type=float, default=1e-3,
                        help="Error rate in sequencing.")
    parser.add_argument('-pc', '--pcr-cycles', type=int, default=0,
                        help="Number of PCR cycles to simulate")
    parser.add_argument('-pe', '--pcr-efficency', type=float, default=1.0,
                        help="PCR amplfication efficency.")
    parser.add_argument('-p', '--processes', type=int, default=1,
                        help="Multiprocessing. Number of processes to run.")
    parser.add_argument('--debug', help="Debug mode", default=False, action="store_true")
    parser.add_argument('-o', "--output", help="Outout file with true/false information.", type=str, default=None)
    parser.add_argument('-of', "--output-format", help="Outout file format. Default: cd-hit", type=str, default="cd-hit")

    return parser.parse_args()


def create_barcode(sequence_options):
    generated_barcode = str()
    for base_options in sequence_options:
        base = random.choice(base_options)
        generated_barcode += base

    logger.debug(f"Barcode created: {generated_barcode}")
    true_barcodes.add(generated_barcode)
    return generated_barcode


def create_barcodes_generator(n, barcode_options, no_prog=False):
    for i in range(n):
        counter['Barcodes start'] += 1
        yield create_barcode(barcode_options)


def amplify(barcode, efficiency=1, error_rate=0.00001):
    """
    simulate amplification on a single barcode efficiency
    """
    def add_base(base, error_rate=0.00001):
        if np.random.random() <= error_rate:
            counter['pcr errors'] += 1
            return np.random.choice(errors_dict[base])
        else:
            return base

    def get_product(barcode, error_rate=0.00001):
        if error_rate == 0:
            return barcode, barcode
        else:
            new_barcode = ''.join(add_base(base, error_rate=error_rate) for base in barcode)
            return barcode, new_barcode

    assert 0 < efficiency <= 1, "PCR efficiency must be between 0 and 1"

    logger.debug(f"Barcode amplification: {barcode}")

    if efficiency == 1:
        return get_product(barcode, error_rate=error_rate)
    elif np.random.random() > (1 - efficiency):
        return get_product(barcode, error_rate=error_rate)
    else:
        return barcode,


def simulate_pcr_cycle(barcodes, efficiency=1, error_rate=0.00001, no_prog=True):
    amplicons = []
    for barcode in tqdm(barcodes, desc="Amplifying", unit="bc", disable=no_prog, leave=False):
        amplicons.extend(amplify(barcode, efficiency=efficiency, error_rate=error_rate))
    return amplicons


def pcr_func(args):
    return simulate_pcr_cycle(*args)


def pcr_cycles(barcodes, efficiency=1, pcr_cycles=0, error_rate=0.00001, nprocs=1):
    if pcr_cycles == 0:
        logger.debug(f"No PCR")
        return collections.Counter(barcodes)

    for cycle in tqdm(range(0, pcr_cycles), desc="Running PCR", unit="cycles", leave=False):
        logger.debug(f"PCR cycle {cycle}")
        if nprocs == 1:
            post_cycle = simulate_pcr_cycle(barcodes, efficiency=efficiency, error_rate=error_rate, no_prog=False)
        else:
            chunks = get_chunks(barcodes, size=10000)
            tasks = [(chunk, efficiency, error_rate) for chunk in chunks]
            post_cycle = []
            with multiprocessing.Pool(nprocs) as pool:
                result = pool.imap_unordered(pcr_func, tasks) # tqdm(pool.imap_unordered(pcr_func, tasks), total=len(barcodes), desc="Amplifying", unit="bc", leave=False)
                for r in result:
                    post_cycle.extend(r)

        barcodes = post_cycle

    return collections.Counter(barcodes)


def add_sequencing_errors(barcode_counts, seq_error_rate=0.01):
    """
    Takes a barcode counter object and adds errors at random to simulate sequencing errors
    :param barcode_counts:
    :param seq_error_rate:
    :return:
    """
    new_list = []
    for barcode in tqdm(barcode_counts.elements(), total=sum(barcode_counts.values()), desc="Sequencing BCs", unit="bc", leave=False):
        new_barcode = str()
        for base in barcode:
            if np.random.random() > seq_error_rate:
                new_barcode += base
            else:
                counter['sequencing errors'] += 1
                new_barcode += np.random.choice(errors_dict[base])

        new_list.append(new_barcode)
    return collections.Counter(new_list)


def main():
    args = get_arguments()
    logging.basicConfig(level=logging.INFO if not args.debug else logging.DEBUG,
                        format="%(levelname)s: %(message)s")

    if args.reads:
        args.pcr_cycles = math.ceil(math.log(args.reads/args.number - 2*args.pcr_efficency, 2*args.pcr_efficency))

    #
    # Header
    #

    print('*' * WIDTH)
    print('SIMULATE BARCODE GENERATION, PCR AND SEQUENCING.')
    print('*'*WIDTH)

    print('Command line options:')

    if args.reads:
        args.pcr_cycles = math.ceil(math.log(args.reads/args.number, 2*args.pcr_efficency))
        print('Note: Calculating PCR cycles based on reads!')

    print('-' * WIDTH)
    arguments = [f"{a}: {v}" for a, v in vars(args).items()]
    print("\n".join(arguments))
    print('-' * WIDTH)

    #
    # Create barcodes
    #
    start = time.time()

    barcode_options = [translate(base) for base in args.sequence]

    logging.info(f"Creating barcodes")
    barcodes = create_barcodes_generator(args.number, barcode_options)

    #
    # Run PCR
    #
    start_pcr = time.time()
    logging.info(f"Running PCR")
    final_barcodes = pcr_cycles(barcodes,
                                efficiency=args.pcr_efficency,
                                pcr_cycles=args.pcr_cycles,
                                error_rate=args.error_rate_pcr,
                                nprocs=args.processes)

    logging.info(f"PCR done, time: {time.time() - start_pcr:.2f} s")
    logging.info(f"PCR errors generated: {counter['pcr errors']:,}")

    #
    # Sequnencing
    #
    start_seq = time.time()

    logging.info(f"Sequencing")

    barcodes_after_seq = add_sequencing_errors(final_barcodes,
                                               seq_error_rate=args.error_rate_seq)

    logging.info(f"Sequencing done, time: {time.time() - start_seq:.2f} s")
    logging.info(f"Sequencing errors generated: {counter['sequencing errors']:,}")


    #
    # Output
    #
    print('-' * WIDTH)
    print('Results')
    print('-' * WIDTH)
    print(f"Number of barcodes before sequencing:      {counter['Barcodes start']:7,}")
    print(f"Number of uniq barcodes before sequencing: {len(true_barcodes):7,}")
    print(f"Number of barcodes after PCR:              {len(final_barcodes):7,}")
    print(f"Number of barcode molecules after PCR:     {sum(final_barcodes.values()):7,}")
    print(f"Number of barcodes after sequencing:       {len(barcodes_after_seq):7,}")
    print(f"Number of barcode reads after sequencing:  {sum(barcodes_after_seq.values()):7,}")
    print('-' * WIDTH)

    logging.info(f"Total run time: {time.time() - start:.2f} s")

    distribution = sorted(collections.Counter(barcodes_after_seq.values()).items())
    print(f"{'Freq':12} {'Reads':12}")
    for reads, freq in distribution:
        print(f"{freq:9} {reads:9}")


    if args.debug:
        for p in sorted(list(collections.Counter(barcodes).items())):
            print(p)
        print()
        for p in sorted(list(final_barcodes.items())):
            print(p)
        print()
        for p in sorted(list(barcodes_after_seq.items())):
            print(p)

    if args.output:
        if args.output_format == 'cd-hit':
            with dnaio.open(args.output, mode='w', fileformat='fasta') as writer:
                for nr, (barcode, count) in enumerate(iter(barcodes_after_seq.items())):
                    is_true = 0
                    if barcode in true_barcodes:
                        is_true = 1

                    record = dnaio.Sequence(f"{is_true}:{nr}:{count}:{barcode}", barcode)
                    writer.write(record)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
