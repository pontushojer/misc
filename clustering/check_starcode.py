"""
Check starcode output for clustering metrics.
"""

import argparse
import pandas as pd
import numpy as np
from collections import Counter, namedtuple
import itertools
import Levenshtein
from tqdm import tqdm
import sys

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("starcode_output", type=str,
                        help="Input starcode output file. Use '-' for stdin. Each line is: "
                             "<MAIN BARCODE> <COUNT> <CLUSTERED BARCODE 1>,<CLUSTERED BARCODE 2>,...")
    parser.add_argument("-t", "--true-barcodes", help="File with true barcodes for simulated data", default=None,
                        type=str)
    parser.add_argument('-o', '--output', help="Output tsv file", default=None)
    parser.add_argument('-s', '--simulated', default=False, action='store_true',
                        help="Based on simulated barcodes,requires file with true_barcodes.")
    parser.add_argument('-d', '--dist-metrics', default=False, action='store_true',
                        help="Include dist metrics for cluster, may take long time.")

    return parser.parse_args()


def main():
    args = get_arguments()

    #dist_func = hamming_distance
    dist_func = Levenshtein.distance

    # If
    barcodes_true =set()
    if args.simulated:
        if not args.true_barcodes:
            sys.exit('Need true barcodes for simulated reads analysis')

        with open(args.true_barcodes, "r") as file:
            barcodes_true = {barcode.strip() for barcode in file.readlines()}

    clusters = list()
    count_categorys = Counter()
    for cluster in tqdm(yield_cluster(args.starcode_output), desc="Reading clusters", unit="clstr"):

        if args.simulated:
            cluster_components_true = cluster.components.intersection(barcodes_true)

            if len(cluster_components_true) == 1:
                cluster.category = "True cluster"
            elif len(cluster_components_true) > 1:
                cluster.category = "Merge"
            else:
                cluster.category = "False cluster"

        count_categorys[cluster.category] += 1

        if args.dist_metrics:
            cluster.get_dist(dist_func)

        clusters.append(cluster)

    print(f"Total clusters: {len(clusters)}")
    print(f"Total true barcodes: {len(barcodes_true)}")
    for cat in count_categorys:
        print(f"Category {cat}: {count_categorys[cat]}")

    df = pd.DataFrame.from_records([cluster.to_dict() for cluster in clusters]).set_index('barcode')

    print(f"Mead cluster read count: {np.mean(df['count']):.2f}")
    print(f"Mead cluster size: {np.mean(df['size']):.2f}")

    if args.dist_metrics:
        external_dists = []
        mean = 0
        for seq1, seq2 in tqdm(itertools.combinations(list(df.index.values), r=2),
                               total=int((len(clusters)*(len(clusters)-1))/2), desc="Reading pairs", unit="pairs"):
            dist = dist_func(seq1, seq2)
            if dist:
                external_dists.append(dist)

            # Check every 1000,000 pairs if change within threshold and break in that case
            if len(external_dists) % 1000_000 == 0:
                current_mean = np.mean(external_dists)

                if abs(mean-current_mean) < 0.001:
                    break
                else:
                    mean = current_mean

        mean_external_dist = np.mean(external_dists)
        mean_internal_dist = np.mean(df.dist)
        print(f'mean_external_dist = {mean_external_dist}')
        print(f'mean_internal_dist = {mean_internal_dist}')

    if args.output:
        df.to_csv(args.output, sep='\t')


def yield_cluster(file):
    if file == '-':
        reader = sys.stdin
    else:
        reader = open(file, "r")

    for line in reader:
        main_barcode, count, barcodes = line.strip().split("\t")

        barcodes = set(barcodes.split(','))

        barcodes.add(main_barcode)

        yield Cluster(main_barcode, count, barcodes)

    reader.close()


def hamming_distance(s1, s2, error_return=None):
    # From https://pythonadventures.wordpress.com/2010/10/19/hamming-distance/
    try:
        assert len(s1) == len(s2)
    except AssertionError:
        return error_return

    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


class Cluster(object):
    def __init__(self, barcode, count, components):
        self.barcode = barcode
        self.count = int(count)
        self.components = components
        self.size = len(self.components)
        self.category = "Unknown"
        self.dist = None

    def get_dist(self, dist_func):
        dists = [dist_func(self.barcode, bc) for bc in self.components if bc is not self.barcode]
        self.dist = np.mean(dists)

    def to_dict(self):
        return {
            "barcode": self.barcode,
            "count": self.count,
            "size": self.size,
            "category": self.category,
            "dist": self.dist
        }


if __name__ == "__main__":
    main()
