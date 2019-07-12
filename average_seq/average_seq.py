import dnaio
import argparse
import sys
import pandas as pd
import matplotlib.pyplot as plt


def get_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("file", help="Input fasta or fastq file", type=str)
    parser.add_argument("-L", "--limit", help="Limit percent to call base", default=0.9, type=float)
    parser.add_argument("-p", "--plot", help="Plot base percents", default=False, action="store_true")
    parser.add_argument("-s", "--subplots", help="Plot as subplots", default=False, action="store_true")
    parser.add_argument("-g", "--gap-char", help="Which charcter to fill the gap with. Default='_'", default="_",
                        type=str)
    parser.add_argument("-o", "--output", help="Output file with the base % for each position.", default=None,
                        type=str)
    return parser.parse_args()


def main():
    args = get_arguments()

    base_list = []

    base_dict = {"A": 0,
                 "T": 0,
                 "C": 0,
                 "G": 0}

    with dnaio.open(args.file, mode="r") as reader:
        for nr, read in enumerate(reader):

            seq = read.sequence

            seq_bases = list(seq)

            create_list = False

            if not base_list:
                create_list = True
            elif len(base_list) != len(seq):
                if len(base_list) > len(seq):
                    seq_bases += ["N"] * (len(base_list) - len(seq))
                else:
                    seq_bases = seq_bases[:len(base_list)]

            for i, base in enumerate(seq_bases):
                if create_list:
                    base_count = base_dict.copy()
                    count_base(base_count, base)
                    base_list.append(base_count)
                else:
                    count_base(base_list[i], base)

    print('Sequence:')

    for entry in base_list:
        total = sum(entry.values())
        for base, count in entry.items():
            if count / total >= args.limit:
                print(base, end="")
                break
        else:
            print(args.gap_char, end="")

    print('')
    percents = {"A": [],
                "T": [],
                "C": [],
                "G": [],
                "Position": []}

    for position, entry in enumerate(base_list, 1):
        total = sum(entry.values())
        percents["Position"].append(position)
        for base, count in entry.items():
            percents[base].append(count / total)

    df = pd.DataFrame(percents, columns=sorted(percents.keys())).set_index("Position")
    print(df)

    if args.output:
        df.to_csv(args.output, sep="\t")

    if args.plot or args.subplots:
        df.plot(kind="line", subplots=args.subplots)
        plt.show()


def count_base(base_dict, base):
    translate = {"A": ["A"],
                 "T": ["T"],
                 "C": ["C"],
                 "G": ["G"],
                 "R": ["A", "G"],
                 "Y": ["C", "T"],
                 "M": ["A", "C"],
                 "K": ["G", "T"],
                 "S": ["C", "G"],
                 "H": ["A", "C", "T"],
                 "B": ["C", "G", "T"],
                 "D": ["A", "G", "T"],
                 "V": ["A", "C", "G"],
                 "N": ["A", "G", "C", "T"]}

    for result_base in translate[base]:
        base_dict[result_base] += 4/len(translate[base])


if __name__ == "__main__":
    main()
