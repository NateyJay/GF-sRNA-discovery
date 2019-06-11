# GF-sRNA-discovery
Clustering and condensing sRNAs to their most expressed variants without use of a genome.

### Motivation
This project was motivated by an urge to do sRNA-seq analysis in species where there is no reference assembly available. Traditionally this is very difficult, as all sRNA loci are known to produce a spectrum of sequence variants, caused by:
1. Loci which do not have specific processing (e.g. heterochromatic siRNAs, long hairpin siRNAs).
2. Imperfect dicing by processing enzymes (DCL proteins), which can produce different sizes and shifted products.
3. Sequencing or transcription errors, making infrequent variants which are functionally the same sRNA.

Reference assemblies overcome this issue by stitching these variants into genetic loci. A genome-free method has to do this by only comparing reads to themselves.

### Our solution
To simplify this complex set of sRNAs, this pipeline condenses variants within an edit_distance to the sequence of the most expressed variant, retaining the abundance levels for tissues. This process is divided into two scripts:
* **01-depth_stacking.py** - Reads in sequence libraries (FASTA), derives abundance for each library for a given read sequence, and then ranks the reads in order of abundance
* **02-clustering_by_kmers.py** - Uses a kmer based approach to cluster reads based on similarity in a greedy manner, working in order of decreasing abundance.

### Usage
A standard usage of these scripts runs them in a linear order, first **depth_stacking** followed by **clustering**.

Program calls:

    # This step is relatively fast and produces intermediate files
    python 01-depth_stacking.py \
      -name [identifier] \        # name for data set
      -lib_files [file_1.fa] [file_2.fa] [file_3.fa] [...] \
      -lib_names [name_1] [name_2] [name_3] [...]
      
    # This step is very time consuming. Expect multi-hour runtimes.
    python 02-clustering_by_kmers.py \
      -name [identifier]\        # name for data set
      -k [integer, default 10]\  # Kmer size for finding similar readsize, larger = faster but less sensitive
      -rpm [float, default 0.5]  # minimum rpm depth for a read to be considered as a centroid
      
Example:

    python 01-depth_stacking.py \
      -name ccm_x_ath \
      -lib_files ./PS_1.fa ./PS_2.fa ./IN_1.fa ./IN_2.fa \
      -lib_names PS1 PS2 IN1 IN2
    
    python 02-clustering_by_kmers.py -name ccm_x_ath
    
### Drawbacks
1. These scripts really only work for 20-24 nt sRNAs and are not optimised nor will function for looking at other size ranges.
2. Modification of the kmer word-size should be done at your own risk. This may cause loss of sensitivity if you raise it above 10. Lowering the word-size will not have any expected advantages for sRNAs.
3. **The biology of this is challenging**. This process makes it easier to analyze these reads in a pinch, but is likely collapsing disparate elements together. For example: miRNA paralogs which have a 1-2 nt difference between them would likely be combined as a single miRNA under this method.
    
### Required Packages
1. python-Levenshtein (0.12.0 tested) [https://github.com/ztane/python-Levenshtein/]
2. tqdm [https://github.com/tqdm/tqdm]
