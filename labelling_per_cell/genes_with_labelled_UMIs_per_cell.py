import pdb
import sys
import os
from subprocess import check_output

umi_stats_file = sys.argv[1]
nMut = int(sys.argv[2])

ofn = "labelled_UMIs_per_cell_" + str(nMut) + "_mutations.csv"

cell_dict = {}

print("estimating file size...")
linecount = int(str(check_output(["wc", "-l", umi_stats_file])).split("'")[1].split(" ")[0])

#pdb.set_trace()

print("processing...")
with open(umi_stats_file) as umi_stats:
    next(umi_stats)
    lines_processed = 0
    for line in umi_stats:
        lines_processed+=1
        if lines_processed % 100000 == 0:
            sys.stdout.write("\r" + str(round(lines_processed/linecount*100, 2)) + " %")
        try:
            lin = line.split("\t")
            CB = lin[1]
            gene = lin[3]
            if not CB in cell_dict:
                cell_dict[CB] = []
                # first dict has labelled gene names, second dict has all gene names
                cell_dict[CB].append([])
                cell_dict[CB].append([])
            if not gene in cell_dict[CB][1]:
                cell_dict[CB][1].append(gene)
            if int(lin[5]) >= nMut:
                if not gene in cell_dict[CB][0]:
                    cell_dict[CB][0].append(gene)
        except KeyError:
            pdb.set_trace()
print("\r100.00%\nwriting statistics...")

with open(ofn, "w") as outf:
    header = "CB\tNrLabelledGenes\tLabelledGenes\tNrTotalGenes\tAllGenes\n"
    outf.write(header)
    for entry in cell_dict:
        line = entry + "\t" + str(len(cell_dict[entry][0])) + "\t" + ",".join(cell_dict[entry][0]) + "\t" + str(len(cell_dict[entry][1])) + "\t" + ",".join(cell_dict[entry][1])+"\n"
        outf.write(line)
