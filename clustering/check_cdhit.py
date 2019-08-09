"""

"""
# cd-hit-454 -i $file -o $file'.clustered' -T $processors -c 0.9 -gap 100  -g 1 -n 3 -M 0

import argparse
import pandas as pd
import numpy as np
from collections import Counter
import itertools
import Levenshtein
from tqdm import tqdm
import sys

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("clstr", help="Input .clstr file. Use '-' for stdin.")
    parser.add_argument('-o', '--output', help="Output tsv file", default=None)
    parser.add_argument('-s', '--simulated', default=False, action='store_true',
                        help="Based on simulated barcodes, first value in name string will tell if true bc or not. "
                             "\nExample)\n"
                             "'>1:12:53:ATCGCTGCATGCTAGCTA' is a true barcode\n"
                             "'>0:12:53:ATCGCTGCATGCTAGCTA' is not")
    parser.add_argument('-d', '--dist-metrics', default=False, action='store_true',
                        help="Include dist metrics for cluster, may take long time.")
    return parser.parse_args()


def main():
    args = get_arguments()

    dist_func = Levenshtein.distance

    clusters = list()
    barcodes_true = 0
    count_categorys = Counter()
    for cluster_id, cluster_list in tqdm(yield_cluster(args.clstr), desc="Reading clusters", unit="clstr"):
        cluster = process_cluster(cluster_list, cluster_id, simulated=args.simulated)
        count_categorys[cluster.category] += 1

        if args.simulated:
            barcodes_true += cluster.true_barcodes

        if args.dist_metrics:
            cluster.get_dist(dist_func)

        clusters.append(cluster)

    print(f"Total clusters: {len(clusters)}")
    print(f"Total true barcodes: {barcodes_true}")

    for cat in sorted(list(count_categorys)):
        print(f"Category {cat}: {count_categorys[cat]}")

    df = pd.DataFrame.from_records([cluster.to_dict() for cluster in clusters]).set_index('barcode')
    df = df.sort_values("count", ascending=False)
    print(df.head(10))

    print(f"Mead cluster read count: {np.mean(df['count']):.2f}")
    print(f"Mead cluster size: {np.mean(df['size']):.2f}")

    if args.dist_metrics:
        external_dists = [dist_func(seq1, seq2) for seq1, seq2 in itertools.combinations(df.seq, r=2)]
        mean_external_dist = np.mean(external_dists)
        mean_internal_dist = np.mean(df.mean_dist)
        print(f'mean_external_dist = {mean_external_dist}')
        print(f'mean_internal_dist = {mean_internal_dist}')

    if args.output:
        df.to_csv(args.output, sep='\t')


def yield_cluster(clstr_file):
    cluster = list()
    cluster_id = int()

    if clstr_file == '-':
        reader = sys.stdin
    else:
        reader = open(clstr_file, "r")

    for line in reader:
        if line.startswith(">"):
            if cluster:
                yield cluster_id, cluster

            cluster = list()
            cluster_id += 1
        else:
            cluster.append(line)

    yield cluster_id, cluster
    reader.close()


def process_cluster(cluster, id_nr, simulated=False):
    """
    :param cluster: List of cluster components.
    :param id_nr: Id of cluster
    :param simulated: boolean
    """
    # Reads cluster file and saves as dict
    cluster_id = int()

    barcodes = set()
    repr_sequence = {"seq": str(), "count": 0}

    c = Counter()
    #read_per_barcode = list()
    for line in cluster:
        barcode = line.strip().split()[2].replace('>', '').replace('.', '')
        if simulated:
            is_true, barcode_nr, barcode_count, barcode_sequence = barcode.split(":")
        else:
            barcode_nr, barcode_count, barcode_sequence = barcode.split(":")
            is_true = 0

        c["Reads"] += int(barcode_count)
        c['True barcodes'] += int(is_true)

        #read_per_barcode.append(int(barcode_count))
        barcodes.add(barcode_sequence)

        # Selects the most numerous sequence as the representative one.
        if repr_sequence['count'] < int(barcode_count):
            repr_sequence['seq'] = barcode_sequence
            repr_sequence['count'] = int(barcode_count)

    cluster = Cluster(repr_sequence['seq'], c['Reads'], barcodes, c['True barcodes'])

    if simulated:
        if c['True barcodes'] == 1:
            cluster.category = "True cluster"
        elif c['True barcodes'] > 1:
            cluster.category = "Merge"
        else:
            cluster.category = "False cluster"

    return cluster


def hamming_distance(s1, s2, error_return=None):
    # From https://pythonadventures.wordpress.com/2010/10/19/hamming-distance/
    try:
        assert len(s1) == len(s2)
    except AssertionError:
        return error_return

    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


class Cluster(object):
    def __init__(self, barcode, count, components, true_barcodes):
        self.barcode = barcode
        self.count = int(count)
        self.components = components
        self.true_barcodes = true_barcodes
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
