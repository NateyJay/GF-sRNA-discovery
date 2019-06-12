import sys
import os
from tqdm import tqdm
import Levenshtein
import argparse
import time
import itertools


parser = argparse.ArgumentParser()
parser.add_argument('-name', nargs='?',required=True)
parser.add_argument('-k', nargs='?', default=10)
parser.add_argument('-rpm', nargs='?', default=0.5)
args = parser.parse_args()
# print args
organismal_suffix = args.name
k_value = int(args.k)
rpm = float(args.rpm)


print ""
print "# Clustering sRNAs by sequence similarity"
print "#",time.ctime(int(time.time()))
print ""
print "            file_id: {0}".format(organismal_suffix)
print "            k_value: {0}".format(k_value)
print "centroid RPM cutoff: {0} rpm".format(rpm)
print "      edit_distance: 2"
print ""

# Kmer clustering algorithm

# This approach uses Kmers to heuristically make edit distance comparisons through sequences in a abundance sorted input file.

# The Kmers approach has the benefit of allowing us to test distances between reads which have a k-word size in common, rather than all reads regardless of similarity. This greatly reduces the number of string distance comparisons compared to what would be needed in an all-by-all method.

# Despite the hueristic approach, this process still takes an enourmous amount of time (2-3 days). Pbar and tqdm provide good estimation of the required time to complete the process. Higher word size (k_value) will greatly increase the speed of this method, but would correlate with increasing false negatives for finding cluster-able reads. That being said, a K of 8, 9 or 10  could be acceptible for sRNAs, while a K of 7 should have no false negatives.



# This function finds a list of words of word length K that make up the a read.
def find_kmers(seq, k_value):

	kmers = []
	longseq = seq + seq
	for i in range(len(seq)):
		kmers.append(longseq[i:i+k_value])


	kmers = list(set(kmers))
	return(kmers)


# reading in all reads in sorted and stacked libraries

dir_path = os.path.dirname(os.path.realpath(__file__))
sorted_loc = '01out-depth_sorted.%s.txt' % (organismal_suffix)

print "Reading in lines..."
# stale_dict = {}
lines = []
with open(sorted_loc, "r") as f:
	lines = f.readlines()

header = lines[0].strip().split()
lines = lines[1:]
# lines = lines[:10000]
total = len(lines)
# sys.exit()


# measures library depth to ascertain the abundance for the RPM cutoff

print "Assessing library depth..."

splines = []
total_abd = 0
for line in lines:
	spline = line.strip().split()

	total_abd += int(spline[-1])
	splines.append(spline)

lines = list(splines)
del splines

print ''
print "found {0:,} unique sequences out of {1:,} total reads".format(len(lines), total_abd)
# print total_abd
# print "rpm 5.0", total_abd / 5000000

rpm_cutoff = total_abd / 1000000 * rpm
print "{0} read per million: {1} reads".format(rpm, rpm_cutoff)
# print "rpm 0.5", total_abd / 500000
print ''


print "Trimming away unique sequences < {0} abd from centroid consideration...".format(rpm_cutoff)
centroid_candidates = []
for line in lines:
	abd = line[-1]
	if int(abd) < rpm_cutoff:
		break
	else:
		centroid_candidates.append(line)


print "{0:,} candidates for cluster centroids".format(len(centroid_candidates))
print ""
# sys.exit()


# here the khash is produced, providing a lookup dictionary for all reads which contain the k-word. 
# this starts by making a dict entry for every possible k-word.

print "Making empty khash..."
khash = {}
DNA = "NATGC"
for output in itertools.product(DNA, repeat=k_value):
	kmer = "".join(output)
	khash[kmer] = []

# finds k-words for all reads and populates khash

times = []
print "Building hash indexes..."
pbar = tqdm(total=total)

inverse_khash = {}

for index, line in enumerate(lines):
	pbar.update(1)
	seq = line[1]
	kmers = find_kmers(seq, k_value)

	for k in kmers:
		khash[k].append(index)

pbar.close()
lookup_lines = list(lines)

# removes khashes which are not found in reads

print "Cleaning-up empty khashes..."
for key in khash.keys():
	entry = khash[key]
	if entry == []:
		del khash[key]




# This step throughs out al reads with an abundance of less than 5 as a cluster. They can still contribute to clusters, but no cluster may have a most expressed sequence of less than 5 reads. Saves this as a list (centroid candidates).


print "Clustering by edit distance..."

structure= {}
c_list= []

tot = len(centroid_candidates)
pbar = tqdm(total=tot)

for centroid_counter, spline in enumerate(centroid_candidates):
	pbar.update(1)

	head_read = spline[0]
	head_seq = spline[1]


	if lookup_lines[centroid_counter] != 0:

		to_remove = find_kmers(head_seq, k_value)
		
		for k in to_remove:
			# print centroid_counter
			# print khash[k]
			khash[k].remove(centroid_counter)


		lookup_lines[centroid_counter] = 0

		c_list.append(centroid_counter)
		structure[centroid_counter] = [centroid_counter]

		kmers = find_kmers(head_seq, k_value)

		found = []

		for k in kmers:
			for i in khash[k]:
				found.append(i)
		found = list(set(found))

		# print "{:,} candidates".format(len(found))
		true_founds = 0
		for i in found:
			if lookup_lines[i] != 0:

				temp_line = lines[i]

				temp_read = temp_line[0]
				temp_seq = temp_line[1]

				distance = Levenshtein.distance(head_seq, temp_seq)

				if distance <= 2:	
					true_founds += 1
					structure[centroid_counter].append(i)
					lookup_lines[i] = 0

					to_remove = find_kmers(temp_seq, k_value)

					# del inverse_khash[i]
					for k in to_remove:
						khash[k].remove(i)
		# print "{:,} true".format(true_founds)
		# print ''
		# if centroid_counter == 5:
		# 	sys.exit()


pbar.close()

# once the structure of the clustering is found, the program writes 2 documents: 
# 1) the list of centroids with abundances
# 2) the list of all reads, with their centroids also indicated

output_loc = '%s/02out-k%s.%s.clustered.txt' % (dir_path, k_value, organismal_suffix)
with open(output_loc, "w") as f:
	print >> f, "type\tcluster\tread\tseq\tlen\t{}".format("\t".join(header[2:]))
centroid_loc = '%s/02out-k%s.%s.centroids.txt' % (dir_path, k_value, organismal_suffix)
with open(centroid_loc, "w") as f:
	print >> f, "cluster\tread\tseq\tlen\t{}\t20\t21\t22\t23\t24".format("\t".join(header[2:]))

print ''
print "{0:,} centroids discovered".format(len(c_list))
print ''
print "Writing to output files..."

tot = len(c_list)
pbar = tqdm(total=tot)
centroid_counter = 0
col_range = range(2, len(header))
for c_i in c_list:
	pbar.update(1)
	centroid_counter += 1

	c_name = "Cl_%s_%s" % (organismal_suffix, centroid_counter)

	reads_in_c = structure[c_i]

	len_abd = {}
	for i in [20,21,22,23,24]:
		len_abd[i] = 0
	lib_abd = {}
	for i in col_range:
		lib_abd[i] = 0


	out_reads = []

	for r_i in reads_in_c:
		full_read = lines[r_i]


		for i in col_range:
			lib_abd[i] += int(full_read[i])

		read_len = len(full_read[1])

		len_abd[read_len] += int(full_read[-1])
		# sys.exit()
		out_reads.append(full_read)

	c_read = lines[c_i]
	c_seq = c_read[1]
	full_c = [c_name] + [c_read[0]] +[c_seq] + [len(c_seq)]

	for i in col_range:
		full_c = full_c + [lib_abd[i]]

	for i in [20,21,22,23,24]:
		full_c = full_c + [len_abd[i]]

	with open(centroid_loc, 'a') as f:

		to_print = "\t".join(map(str, full_c))
		print >> f, to_print

	# print len(out_reads)
	with open(output_loc, 'a') as f:
		to_print = ["C"] + full_c[:-5]
		to_print = "\t".join(map(str, to_print))
		print >> f, to_print

		for read in out_reads:
			length = len(read[1])
			to_print = ["H", c_name] + read
			to_print.insert(4,length)
			# print to_print
			to_print = "\t".join(map(str, to_print))
			print >> f, to_print


pbar.close()




	# sys.exit()







