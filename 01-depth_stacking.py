import sys
import operator
import argparse
import os


# reading in arguments
parser = argparse.ArgumentParser()
parser.add_argument('-name', nargs='?',required=True)
parser.add_argument('-lib_files', nargs='+', required=True)
parser.add_argument('-lib_names', nargs='+', required=True)
args = parser.parse_args()

files = args.lib_files
libs = args.lib_names
org = args.name





dir_path = os.path.dirname(os.path.realpath(__file__))

# this is the key for counting reads specific to libraries
def make_key(seq, lib):
	return((seq, lib))
	# return("%(seq)s-%(lib)s" % locals())

print "Processing:", org
print "Parsing to dict..."


line_count = 0
perc_count = 0
seqs = []

acceptable_sizes = set([20,21,22,23,24])


# reads through files, maintaining a count of each read for each library
read_dict = {}

for i, file in enumerate(files):
	lib = libs[i]
	print "  Reading:", lib, file

	new_key_count = 0

	with open(file, "r") as f:
		for line in f:

			if line[0] != ">":
				seq = line.strip()

				if len(seq) in acceptable_sizes:
				
					key = make_key(seq, i)

					try:
						read_dict[key] += 1

					except KeyError:	
						new_key_count += 1
						# for i in libs:
						# 	new_key = make_key(seq, i)
						# 	read_dict[new_key] = 0
						seqs.append(seq)
						read_dict[key] = 1
	print "  {:,} new keys added".format(new_key_count)


# produces a list of reads and total abundance in all libraries, which is sorted in descending order

print "Sorting by depth..."
reads = []
# keys = keynames(read_dict)
for seq in seqs:
	entry = [seq]
	abd = 0
	for i, lib in enumerate(libs):
		
		key = make_key(seq, i)

		try:
			lib_abd = read_dict[key]
			del read_dict[key]
		except KeyError:
			lib_abd = 0

		entry.append(lib_abd)
		abd += lib_abd
	entry.append(abd)
	reads.append(entry)


reads.sort(key=operator.itemgetter(len(reads[0])-1))


# initializes and writes the header line of the document

stacked_loc = "{}/01out-depth_sorted.{}.txt".format(dir_path,org)
with open(stacked_loc, "w") as file:
	header = ["read","seq"]
	for lib in libs:
		header.append(lib)
	header.append("tot_abd")

	print >> file, "\t".join(header)


# processes through lines, retaining reads 20-24 in length and printing their abundances to a file
total = len(reads)
print "  {} lines to process".format(total)
perc_count = 0
line_count = 0
print "Writing to: %s" % (stacked_loc)

with open(stacked_loc, "a") as stacked_file:


	for read in reads[::-1]:
		if line_count % (total / 10) == 0:
			print perc_count
			perc_count += 10
		line_count += 1

		seqlen = len(read[0])
		if seqlen >= 20 and seqlen <=24:


			read_name = "R_%s" % (line_count)

			to_print = [read_name] + read
			to_print = "\t".join(map(str, to_print))

			print >> stacked_file, to_print








