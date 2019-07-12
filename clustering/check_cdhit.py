
import argparse
import pandas as pd
from collections import Counter

def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("clstr", help="Input .clstr file")
    parser.add_argument('-o', '--output', help="Output tsv file", default=None)

    return parser.parse_args()


def main():
    args = get_arguments()

    clusters = list()
    with open(args.clstr, "r") as clstr_file:
        for cluster_id, cluster_list in yield_cluster(clstr_file):
            clusters.append(process_cluster(cluster_list, cluster_id))

    df = pd.DataFrame(clusters).set_index('id')
    df = df.sort_values(by=["sequences", "reads"], ascending=False)
    counts_seq = Counter(df['sequences'])
    counts_reads = Counter(df['reads'])
    counts_seq = sorted(list(counts_seq.items()))
    counts_reads = sorted(list(counts_reads.items()))

    counts_size = df[["sequences","reads"]].groupby("sequences").sum()
    print(counts_size)

    print(f"Cluster size (sequences), Frequency")
    for bcs, frequency in counts_seq:
        print(f"{bcs:12}, {frequency:9}")

    print(f"Cluster size (reads), Frequency")
    for reads, frequency in counts_reads:
        print(f"{reads:12}, {frequency:9}")


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


def process_cluster(cluster, id_nr):
    """
    :param cluster: List of cluster components.
    :param id_nr: Id of cluster
    """

    # Reads cluster file and saves as dict
    cluster_dict = dict()
    cluster_id = int()

    clusters = []
    c = Counter()
    read_per_barcode = list()
    for line in cluster:
        barcode = line.strip().split()[2].replace('>', '').replace('.', '')
        barcode_nr, barcode_count, barcode_sequence = barcode.split(":")
        c["Unique barcodes"] += 1
        c["Reads"] += int(barcode_count)
        read_per_barcode.append(int(barcode_count))

    return {'id': id_nr,
            'sequences': c["Unique barcodes"],
            'reads': c["Reads"],
            'read_per_bc': sorted(read_per_barcode, reverse=True)}


if __name__ == "__main__":
    main()
