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

## it can now also be used with and without a SNP file

## it now has variable tag format (CR or STAR/Rsubread style)


if len(sys.argv) == 6:
    ifn = sys.argv[1]
    qual_thr = int(sys.argv[2])
    snp_file = sys.argv[3]
    cutoff_rate = int(sys.argv[4])/100
    tag_format = sys.argv[5]
else:
    print("Please provide me with an input .bam file to find mutations in, a minimum quality score for said mutations and a file containing known SNP positions and a cutoff rate for SNPs, if you want to filter out some positions. Otherwise type N and 101. Additionally please let me know if the file uses cellranger [CR] or STAR/Rsubread [ST] gene tags. CB/UB tags are expected for cell barcode and UMI.")
    sys.exit()

if tag_format == "CR" or tag_format == "ST":
    if tag_format == "CR":
        gene_tag = "GN"
        ensmbl_tag = "GX"
    else:
        gene_tag = "XT"
        ensmbl_tag = "XT"
else:
    print("Please indicate the tag format for gene tags - CR for cellranger (GN/GX - tags) or ST for STAR/Rsubread (XT - tags)")
    sys.exit()

## read in SNP file and create a set
snp_set = []
if snp_file != "N":
    with open(snp_file, "r") as snps:
        next(snps)
        for line in snps:
            lin = line.rstrip().split("\t")
            #allowing only those positions in that have a TC rate over the desired cutoff rate
            if int(lin[2]) > 0:
                if int(lin[2]) / int(lin[1]) >= cutoff_rate:
                    snp_set.append(lin[0])
else:
    snp_set = ["0"]
                        #pdb.set_trace()
    #to speed up membership tests
snp_set = set(snp_set)


samfile = pysam.AlignmentFile(ifn, "rb")
#ofn_coverage = ifn[:-4] + "_coverage_per_gene_TC-Q" + str(qual_thr) + "_SNPs.tsv"
ofn_coverage_T_pos = ifn[:-4] + "_coverage_per_T-position_TC-Q" + str(qual_thr) + "_SNPs.tsv"
ofn_meanMedian_per_gene = ifn[:-4] + "_mean_median_cov_TC_per_gene-Q" + str(qual_thr) + "_SNPs.tsv"
#ofn_genes_positions = ifn[:-4] + "_genes_positions.tsv"

gene_T_occurence_dict = {}
# structure: gene : dict{position : [countT, countTC]}

coverage_per_position_dict = {}
# structure : position[coverage, RefBase, numTC]

coverage_per_T_position_dict = {}
# structure : position[coverage, numTC]

coverage_per_gene_dict = {}
# structure : gene : position[countT, numTC]

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
    if (flag == 0 or flag == 16) and read.has_tag(ensmbl_tag) and not re.search("[ID]", read.cigarstring):
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
        gene = read.get_tag(ensmbl_tag)
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
                    if not gene in coverage_per_gene_dict:
                        coverage_per_gene_dict[gene] = {}
                        coverage_per_gene_dict[gene][glob_position] = [1,0]
                    else:
                        if not glob_position in coverage_per_gene_dict[gene]:
                            coverage_per_gene_dict[gene][glob_position] = [1,0]
                        else:
                            coverage_per_gene_dict[gene][glob_position][0]+=1

            except IndexError:
                pass
                #pdb.set_trace()
        MTtag = read.get_tag("MT").split("_")
        if not MTtag[0] == "0-XX":
            for mutation in MTtag:
                mutation = mutation.split("-")
                #pdb.set_trace()
                try:
                    if mutation[1] == label and int(mutation[2]) >= qual_thr and not mutation[0] in snp_set:
                        coverage_per_T_position_dict[mutation[0]][1]+=1
                        coverage_per_gene_dict[gene][mutation[0]][1]+=1
                except IndexError:
                    pass
                    #pdb.set_trace()
#pdb.set_trace()
sys.stdout.write("\nwriting coverage per T position and TtoC mutations...\n")

f2 = open(ofn_coverage_T_pos, "w")
f2.write("position\tcoverage\tnumTC\n")
for entry in coverage_per_T_position_dict:
    f2.write(entry + "\t" + str(coverage_per_T_position_dict[entry][0]) + "\t" + str(coverage_per_T_position_dict[entry][1]) + "\n")
f2.close()

sys.stdout.write("\nwriting mean / median coverage and mean / median TC per gene...\n")
f3 = open(ofn_meanMedian_per_gene, "w")
f3.write("gene\tMeanTcoverage\tmedianTcoverage\tmeanTC\tmedianTC\n")
for gene in coverage_per_gene_dict:
    coverage_vector = []
    TC_vector = []
    for position in coverage_per_gene_dict[gene]:
        coverage_vector.append(coverage_per_gene_dict[gene][position][0])
        TC_vector.append(coverage_per_gene_dict[gene][position][1])
    #pdb.set_trace()
    cov_mean = np.mean(coverage_vector)
    cov_median = np.median(coverage_vector)
    TC_mean = np.mean(TC_vector)
    TC_median = np.median(TC_vector)
    f3.write(gene + "\t" + str(cov_mean) + "\t" + str(cov_median) + "\t" + str(TC_mean) + "\t" + str(TC_median) + "\n")
f3.close()
