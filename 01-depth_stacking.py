import sys
import operator
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-name', nargs='?',required=True)
parser.add_argument('-lib_files', nargs='+', required=True)
parser.add_argument('-lib_names', nargs='+', required=True)
args = parser.parse_args()

files = args.lib_files
libs = args.lib_names
org = args.name

# print INs
# print PSs
# print org
# sys.exit()



dir_path = os.path.dirname(os.path.realpath(__file__))


def make_key(seq, lib):
	return("%(seq)s-%(lib)s" % locals())

print "Processing:", org
print "Parsing to dict..."


line_count = 0
perc_count = 0
seqs = []

# libs = ["IN1","IN2","IN3","PS1","PS2","PS3"]
read_dict = {}

for i, file in enumerate(files):
	lib = libs[i]
	print "  Reading:", lib, file

	with open(file, "r") as f:
		for line in f:

			if line[0] != ">":
				seq = line.strip()
				
				key = make_key(seq, lib)

				try:
					read_dict[key] += 1

				except KeyError:	
					for i in libs:
						new_key = make_key(seq, i)
						read_dict[new_key] = 0
					seqs.append(seq)
					read_dict[key] += 1


print "Sorting by depth..."
reads = []
# keys = keynames(read_dict)
for seq in seqs:
	entry = [seq]
	abd = 0
	for lib in libs:
		key = make_key(seq, lib)
		lib_abd = read_dict[key]

		entry.append(lib_abd)
		abd += lib_abd
	entry.append(abd)
	reads.append(entry)


reads.sort(key=operator.itemgetter(len(reads[0])-1))


stacked_loc = "{}/01out-depth_sorted.{}.txt".format(dir_path,org)
with open(stacked_loc, "w") as file:
	header = ["read","seq"]
	for lib in libs:
		header.append(lib)
	header.append("tot_abd")

	print >> file, "\t".join(header)


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








