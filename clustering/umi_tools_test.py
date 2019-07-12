import dnaio
import argparse
from umi_tools import UMIClusterer
import pandas as pd
from collections import Counter
import sys
import time

def get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("fasta", help="Input fasta file")

    return parser.parse_args()


def main():
    args = get_arguments()

    bc_dict = dict()
    with dnaio.open(args.fasta, mode="r") as file:
        for read in file:
            bc_id, bc_count, bc_seq = read.name.strip('>').split(':')
            bc_dict[bc_seq] = int(bc_count)

    # Based on https://umi-tools.readthedocs.io/en/latest/API.html
    clusterer = UMIClusterer(cluster_method='directional')

    start = time.time()
    clustered_bcs = clusterer(bc_dict, threshold=1)
    end = time.time()

    cluster_lens = [len(c) for c in clustered_bcs]
    count = Counter(cluster_lens)
    count = sorted(list(count.items()))
    print(f"Cluster size, Frequency")
    for bcs, frequency in count:
        print(f"{bcs:12}, {frequency:9}")

    print(f'Time to run: {end-start} s')
    print(f'Length data: {len(bc_dict)}')

if __name__ == "__main__":
    main()


cd-hit-454 -i  -o $file'.clustered' -T $processors -c 0.9 -gap 100 -g 1 -n 3 -M 0 > $path"/cdhit.log