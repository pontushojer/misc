"""
Silmulate randomly generated barcodes.
Based on UMI-Tools pipeline:
https://github.com/CGATOxford/UMI-tools_pipelines/blob/master/notebooks/Simulating_umi_deduping.ipynb
"""

import argparse
import pandas as pd
import sys
import random
import collections
import numpy as np
import logging
import time
from tqdm import tqdm
import math
import multiprocessing
from queue import Queue

logger = logging.getLogger(__name__)

counter = collections.Counter()

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
    return parser.parse_args()


def create_barcode(sequence_options):
    generated_barcode = str()
    for base_options in sequence_options:
        base = random.choice(base_options)
        generated_barcode += base

    logger.debug(f"Barcode created: {generated_barcode}")
    return generated_barcode


def create_barcodes(n, barcode_options, no_prog=False):
    barcodes = list()
    for i in tqdm(range(n), total = n, desc="Create BCs", unit= "bc", leave=False, disable=no_prog):
        counter['Barcodes start'] += 1
        barcodes.append(create_barcode(barcode_options))
    return barcodes


def create_barcodes_generator(n, barcode_options, no_prog=False):
    for i in tqdm(range(n), total = n, desc="Create BCs", unit= "bc", leave=False, disable=no_prog):
        counter['Barcodes start'] += 1
        yield create_barcode(barcode_options)


def amplify(barcode, efficiency=1, error_rate=0.00001):
    """
    simulate amplification on a single barcode efficiency
    """
    assert 0 < efficiency <= 1, "PCR efficiency must be between 0 and 1"
    logger.debug(f"Barcode amplification: {barcode}")

    if efficiency == 1:
        logger.debug(f"Amplified.")
        if error_rate == 0:
            return barcode, barcode
        else:
            new_barcode = str()
            for base in barcode:
                if np.random.random() <= error_rate:
                    new_barcode += np.random.choice(errors_dict[base])
                else:
                    new_barcode += base
            return barcode, new_barcode
    else:
        if np.random.random() > (1 - efficiency):
            logger.debug(f"Amplified.")
            if error_rate == 0:
                return barcode, barcode
            else:
                new_barcode = str()
                for base in barcode:
                    if np.random.random() <= error_rate:
                        new_barcode += np.random.choice(errors_dict[base])
                    else:
                        new_barcode += base
                return barcode, new_barcode
        else:
            logger.debug(f"Not amplified.")
            return barcode,


def simulate_pcr_cycle(barcodes, efficiency=1, error_rate=0.00001):
    new_list = []
    for barcode in tqdm(barcodes, desc="Amplifing BCs",unit= "bc", leave=False):
        new_list.extend(amplify(barcode, efficiency=efficiency, error_rate=error_rate))
    return new_list


def pcr_cycles(barcodes, efficiency=1, pcr_cycles=0, error_rate=0.00001):
    if pcr_cycles == 0:
        logger.debug(f"No PCR")
        return collections.Counter(barcodes)

    for cycle in tqdm(range(0, pcr_cycles), desc="Running PCR", unit="cycles", leave=False):
        logger.debug(f"PCR cycle {cycle}")
        post_cycle = simulate_pcr_cycle(barcodes, efficiency=efficiency, error_rate=error_rate)

        barcodes = post_cycle
    return collections.Counter(barcodes)


def add_sequencing_errors(barcode_counter, seq_error_rate=0.01):
    """
    Takes a barcode counter object and adds errors at random to simulate sequencing errors
    :param barcode_counter:
    :param seq_error_rate:
    :return:
    """
    new_list = []
    for barcode in tqdm(barcode_counter.elements(), total=sum(barcode_counter.values()), desc="Sequencing BCs",unit= "bc", leave=False):
        new_barcode = str()
        for base in barcode:
            if np.random.random() > seq_error_rate:
                new_barcode += base
            else:
                counter['sequencing errors'] += 1
                new_barcode += np.random.choice(errors_dict[base])

        new_list.append(new_barcode)
    return collections.Counter(new_list)


def func(args):
    return create_barcodes(*args)

def process_create(nprocs, number, barcode_options):

    chunks = [int(math.ceil(number / nprocs)) for i in range(nprocs - 1)]
    chunks = chunks + [number - sum(chunks)]

    tasks = [(chunk, barcode_options, True) for chunk in chunks]
    output = []
    with multiprocessing.Pool(nprocs) as pool:
        result = pool.imap_unordered(func, tasks)
        for r in result:
            output.extend(r)

    return output

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
        args.pcr_cycles = math.ceil(math.log(args.reads/args.number - 2*args.pcr_efficency, 2*args.pcr_efficency))
        print('Note: Calculating PCR cycles based on reads!')

    print('-' * WIDTH)
    arguments = [f"{a}: {v}" for a, v in vars(args).items()]
    print("\n".join(arguments))
    print('-' * WIDTH)

    #
    # Create barcodes
    #
    barcode_options = [translate(base) for base in args.sequence]

    logging.info(f"Creating barcodes")
    start_bc_create = time.time()

    barcodes = create_barcodes(args.number, barcode_options)

    logging.info(f"Barcodes created, time: {time.time() - start_bc_create:.2f} s")

    #barcodes = process_create(args.processes, args.number, barcode_options)
    # barcodes = multiprocess_create(args. processes, args.number, args.sequence)

    #
    # Run PCR
    #

    logging.info(f"Running PCR")
    start_pcr = time.time()

    final_barcodes = pcr_cycles(barcodes,
                                efficiency=args.pcr_efficency,
                                pcr_cycles=args.pcr_cycles,
                                error_rate=args.error_rate_pcr)

    logging.info(f"PCR done, time: {time.time() - start_pcr:.2f} s")

    #
    # Sequnencing
    #

    logging.info(f"Sequencing")
    start_seq = time.time()

    barcodes_after_seq = add_sequencing_errors(final_barcodes,
                                               seq_error_rate=args.error_rate_seq)

    logging.info(f"Sequencing done, time: {time.time() - start_seq:.2f} s")
    logging.info(f"Sequencing errors generated: {counter['sequencing errors']}")


    #
    # Output
    #
    print('-' * WIDTH)
    print('Results')
    print('-'*WIDTH)
    print(f"Number of barcodes before sequencing:     {len(barcodes):7,}")
    print(f"Number of barcodes after PCR:             {len(final_barcodes):7,}")
    print(f"Number of barcode molecules after PCR:    {sum(final_barcodes.values()):7,}")
    print(f"Number of barcodes after sequencing:      {len(barcodes_after_seq):7,}")
    print(f"Number of barcode reads after sequencing: {sum(barcodes_after_seq.values()):7,}")
    print('-' * WIDTH)

    logging.info(f"Total run time: {time.time() - start_bc_create:.2f} s")

    if args.debug:
        for p in sorted(list(collections.Counter(barcodes).items())):
            print(p)
        print()
        for p in sorted(list(final_barcodes.items())):
            print(p)
        print()
        for p in sorted(list(barcodes_after_seq.items())):
            print(p)

if __name__ == "__main__":
    multiprocessing.freeze_support()
    main()
