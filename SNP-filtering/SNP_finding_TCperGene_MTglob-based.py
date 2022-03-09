import pysam
import re
import pdb
import sys
import os
from subprocess import check_output
import statistics
import re
import numpy as np

## this makes a list of all positions of T per genes (and their number of occurence)
## it also lists all positions of TtoC and their number of occurence
## it then tabluates average and median coverage for each gene
## it also calculates what fraction of positions is affected by TtoC mutations
## this one is SNP-aware

## this converts all lowercase letters in the reference to uppercase, not sure if this is a good idea


if len(sys.argv) == 4:
    ifn = sys.argv[1]
    qual_thr = int(sys.argv[2])
    snp_file = sys.argv[3]
else:
    print("Please provide me with an input .bam file to find mutations in, a minimum quality score for said mutations and a file containing known SNP positions.")
    sys.exit()

## read in SNP file and create a set
snp_set = []
with open(snp_file, "r") as snps:
    next(snps)
    for line in snps:
        snp_set.append(line.rstrip().split(",")[1].replace('"', ''))
    #to speed up membership tests
    snp_set = set(snp_set)


samfile = pysam.AlignmentFile(ifn, "rb")
#ofn_coverage = ifn[:-4] + "_coverage_per_gene_TC-Q" + str(qual_thr) + "_SNPs.tsv"
ofn_coverage_T_pos = ifn[:-4] + "_coverage_per_T_position_TC-Q" + str(qual_thr) + "_SNPs.tsv"
#ofn_genes_positions = ifn[:-4] + "_genes_positions.tsv"

gene_T_occurence_dict = {}
# structure: gene : dict{position : [countT, countTC]}

coverage_per_position_dict = {}
# structure : position[coverage, RefBase, numTC]

coverage_per_T_position_dict = {}
# structure : position[coverage, numTC]

# for the progress counter
sys.stdout.write("estimating file size...\n")
linecount = int(str(check_output(["samtools", "view", "-c", ifn]).rstrip()).split("'")[1])
reads_processed = 0
sys.stdout.write("processing...\n")

reads_processed = 0


for read in samfile.fetch():
    reads_processed+=1
    flag = read.flag
    # only primary alignments allowed, for now we disinvite InDel-d reads
    # only things that map to a gene are allowed to join
    if (flag == 0 or flag == 16) and read.has_tag("XT") and not re.search("[ID]", read.cigarstring):
    # and not re.search("[ID]", read.cigarstring):
        if reads_processed % 10000 == 0:
            sys.stdout.write("\r" + str(round(reads_processed/linecount*100, 2)) + " %")
        #setting the base we are looking for, based on if the read is forward or reverse
        refBase = "T"
        label = "TC"
        labelBase = "C"
        if flag == 16:
            refBase = "A"
            label = "AG"
            labelBase = "G"
        contig = read.reference_name
        gene = read.get_tag("XT")
        read_sequence = read.query_alignment_sequence
        reference = read.get_reference_sequence().upper()
        seq_quality = read.query_alignment_qualities
        reference_positions = read.get_reference_positions()
        ## let's find all the T
        T_positions = np.where(np.array(list(reference)) == refBase)[0]
        for T_position in T_positions:
            try:
                if seq_quality[T_position] >= qual_thr:
                    glob_position = contig + ":" + str(reference_positions[T_position])
                    if not glob_position in coverage_per_T_position_dict:
                        coverage_per_T_position_dict[glob_position] = [1,0]
                    else:
                        coverage_per_T_position_dict[glob_position][0]+=1
            except IndexError:
                pdb.set_trace()
        MTtag = read.get_tag("MT").split("_")
        if not MTtag[0] == "0-XX":
            for mutation in MTtag:
                mutation = mutation.split("-")
                #pdb.set_trace()
                try:
                    if mutation[1] == label and int(mutation[2]) >= qual_thr and not mutation[0] in snp_set:
                        coverage_per_T_position_dict[mutation[0]][1]+=1
                except IndexError:
                    pdb.set_trace()

sys.stdout.write("\nwriting coverage per T position and TtoC mutations...\n")

f2 = open(ofn_coverage_T_pos, "w")
f2.write("position\tcoverage\tnumTC\n")
for entry in coverage_per_T_position_dict:
    f2.write(entry + "\t" + str(coverage_per_T_position_dict[entry][0]) + "\t" + str(coverage_per_T_position_dict[entry][1]) + "\n")
f2.close()
