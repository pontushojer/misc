"""

"""
# cd-hit-454 -i $file -o $file'.clustered' -T $processors -c 0.9 -gap 100  -g 1 -n 3 -M 0

import argparse
import pandas as pd
import numpy as np
from collections import Counter


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("clstr", help="Input .clstr file")
    parser.add_argument('-o', '--output', help="Output tsv file", default=None)
    parser.add_argument('-s', '--simulated', default=False, action='store_true',
                        help="Based on simulated barcodes, first value in name string will tell if true bc or not. "
                             "\nExample)\n"
                             "'>1:12:53:ATCGCTGCATGCTAGCTA' is a true barcode\n"
                             "'>0:12:53:ATCGCTGCATGCTAGCTA' is not")
    return parser.parse_args()


def main():
    args = get_arguments()

    clusters = list()
    with open(args.clstr, "r") as clstr_file:
        for cluster_id, cluster_list in yield_cluster(clstr_file):
            clusters.append(process_cluster(cluster_list, cluster_id, simulated=args.simulated))

    df = pd.DataFrame(clusters).set_index('id')
    df = df.sort_values(by=["sequences", "reads"], ascending=False)
    counts_seq = Counter(df['sequences'])
    counts_reads = Counter(df['reads'])
    counts_seq = sorted(list(counts_seq.items()))
    counts_reads = sorted(list(counts_reads.items()))

    counts_size = df[["sequences", "reads"]].groupby("sequences").sum()
    print(counts_size)

    # print(f"Cluster size (sequences), Frequency")
    # for bcs, frequency in counts_seq:
    #     print(f"{bcs:12}, {frequency:9}")
    #
    # print(f"Cluster size (reads), Frequency")
    # for reads, frequency in counts_reads:
    #     print(f"{reads:12}, {frequency:9}")

    if args.simulated:
        df_new = df[['true_bc', 'sequences', 'reads']]
        df_new = df_new.groupby('true_bc').aggregate({"reads": [np.sum, np.mean],
                                                      "sequences": [np.sum, np.mean, np.size]})
        print(df_new.head())
        counts_true = Counter(df['true_bc'])
        counts_true = sorted(list(counts_true.items()))
        print(f"True bracodes, Frequency")
        for bcs, frequency in counts_true:
            print(f"{bcs:12}, {frequency:9}")

    if args.output:
        df.to_csv(args.output, sep='\t')


def yield_cluster(clstr_file):

    cluster = list()
    cluster_id = int()
    for line in clstr_file:
        if line.startswith(">"):
            if cluster:
                yield cluster_id, cluster

            cluster = list()
            cluster_id += 1
        else:
            cluster.append(line)

    yield cluster_id, cluster


def process_cluster(cluster, id_nr, simulated=False):
    """
    :param cluster: List of cluster components.
    :param id_nr: Id of cluster
    :param simulated: boolean
    """
    # Reads cluster file and saves as dict
    cluster_dict = dict()
    cluster_id = int()

    clusters = []
    c = Counter()
    read_per_barcode = list()
    for line in cluster:
        barcode = line.strip().split()[2].replace('>', '').replace('.', '')
        if simulated:
            is_true, barcode_nr, barcode_count, barcode_sequence = barcode.split(":")
        else:
            barcode_nr, barcode_count, barcode_sequence = barcode.split(":")
            is_true = 0
        c["Unique barcodes"] += 1
        c["Reads"] += int(barcode_count)
        c['True barcodes'] += int(is_true)
        read_per_barcode.append(int(barcode_count))

    return {'id': id_nr,
            'sequences': c["Unique barcodes"],
            'reads': c["Reads"],
            'true_bc': c['True barcodes'],
            'read_per_bc': sorted(read_per_barcode, reverse=True)}


if __name__ == "__main__":
    main()


