import sys
import csv
import itertools


def hamming_distance(s1, s2):
	# From https://pythonadventures.wordpress.com/2010/10/19/hamming-distance/
	assert len(s1) == len(s2)
	return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def main():
	input = sys.argv[1]

	with open(input,"r") as tsvfile:
		tsvreader = csv.reader(tsvfile, delimiter="\t")
		bcs = []
		for name, sequence in tsvreader:
			bcs.append({"Name": name, "Seq": sequence})

	print (bcs) 
	
	for bc1, bc2 in itertools.combinations(bcs, 2):
		dist = hamming_distance(bc1["Seq"], bc2["Seq"])
		
		print(f"{bc1['Name']} and {bc2['Name']} at dist = {dist}")


if __name__ == "__main__":
	main()
