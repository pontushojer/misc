import dnaio
import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Input fasta or fastq file", type=str)
    parser.add_argument("-r", "--reference", help="Reference sequence to compare to", default=False, type=str)
    parser.add_argument("-s", "--start", help="Start position for comparison", default=0, type=int)
    return parser.parse_args()


def main():
    args = get_arguments()

    end = args.start + len(args.reference)

    args.reference = args.reference.upper()
    l_ref = len(args.reference)
    # N = true sequence
    n_total = 0
    n_position = [0] * l_ref

    # E = error sequence
    e_total = 0
    e_position = [0] * l_ref

    with dnaio.open(args.file, mode="r") as reader:
        for read in tqdm(reader):
            comp_seq = read.sequence[args.start:end]
            # Check identical
            if comp_seq == args.reference:
                n_total += 1
                n_position = list(map(sum, zip(n_position,[1] * l_ref)))
            else:
                e_total += 1
                for i, (comp_base, ref_base) in enumerate(zip(comp_seq, args.reference)):
                    if comp_base == ref_base:
                        n_position[i] += 1
                    else:
                        e_position[i] += 1
    print(n_position)
    print(e_position)


    error_rate_total = e_total / (e_total + n_total)

    t_position = [sum(tup) for tup in zip(n_position, e_position)]
    print(t_position)

    error_rate_position = [e / t for e, t in  zip(e_position, t_position)]

    print(f"Total error rate:       {error_rate_total:6.3%}")
    print(f"Mean error rate:        {sum(error_rate_position)/l_ref:6.3%}")
    print("-"*30)
    for i, er_pos in enumerate(error_rate_position, args.start + 1):
        print(f"Error rate position {i:2}: {er_pos:6.3%} ")



if __name__ == "__main__":
    main()
