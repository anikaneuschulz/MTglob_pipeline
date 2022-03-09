import pysam
import pdb
import sys
import subprocess
import os
import numpy as np
import time

if len(sys.argv) == 3:
    ifn = sys.argv[1]
    subfolder = sys.argv[2]
else:
    print("Please feed me an indexed .bam file to make a count matrix from and the name of the subfolder to place it in.")
    sys.exit()


bamfile = pysam.AlignmentFile(ifn, "rb")

CB_list = np.array([])
GN_list = np.array([])
nonempty_cells = 0

count_dict = {}
# structure : CB {GN[UB]}

outfolder = "/".join(os.path.abspath(ifn).split("/")[:-1]) + "/" + subfolder
ofn = outfolder + "/matrix.mtx"
barcodes = outfolder + "/barcodes.tsv"
features = outfolder + "/features.tsv"


sys.stdout.write("estimating file size...\n")
linecount = int(str(subprocess.check_output(["samtools", "view", "-c", ifn]).rstrip()).split("'")[1])
reads_processed = 0
current_time = time.strftime("%H:%M:%S", time.localtime())
sys.stdout.write("[" + current_time + "]")
sys.stdout.write("processing...\n")

for read in bamfile.fetch():
    reads_processed+=1
    if reads_processed % 10000 == 0:
        sys.stdout.write("\r" + str(round(reads_processed/linecount*100,2)) + "% of reads processed.")
    flag = read.flag
    # only primary alignments allowed
    if flag == 0 or flag == 16:
        if read.has_tag("CB") and read.has_tag("UB") and read.has_tag("GN"):
            CB = read.get_tag("CB")
            UB = read.get_tag("UB")
            GN = read.get_tag("GN")
            #pdb.set_trace()
            if not GN in GN_list:
                GN_list = np.concatenate((GN_list, np.array([GN])))
                GN_index = np.where(GN_list==GN)[0][0]+1
            else:
                GN_index = np.where(GN_list==GN)[0][0]+1
            if not CB in CB_list:
                CB_list = np.concatenate((CB_list, np.array([CB])))
                CB_index = np.where(CB_list==CB)[0][0]+1
            else:
                CB_index = np.where(CB_list==CB)[0][0]+1
            if not CB_index in count_dict:
                count_dict[CB_index] = {}
                count_dict[CB_index][GN_index] = []
                count_dict[CB_index][GN_index].append(UB)
                nonempty_cells+=1
            else:
                if not GN_index in count_dict[CB_index]:
                    count_dict[CB_index][GN_index] = []
                    count_dict[CB_index][GN_index].append(UB)
                    nonempty_cells+=1
                else:
                    if not UB in count_dict[CB_index][GN_index]:
                        count_dict[CB_index][GN_index].append(UB)
                        nonempty_cells+=1

sys.stdout.write("\nwriting...\n")
current_time = time.strftime("%H:%M:%S", time.localtime())
sys.stdout.write("[" + current_time + "]")

command = "mkdir " + outfolder
os.system(command)

with open(ofn, "w") as o:
    o.write("%%MatrixMarket matrix coordinate integer general\n%\n")
    line = str(len(GN_list)) + " " + str(len(CB_list)) + " " + str(nonempty_cells) + "\n"
    o.write(line)
    for CB in count_dict:
        for GN in count_dict[CB]:
            line = str(GN) + " " + str(CB) + " " + str(len(count_dict[CB][GN])) + "\n"
            o.write(line)

with open(features, "w") as f:
    for gene in GN_list:
        f.write(gene + "\n")

with open(barcodes, "w") as f:
    for CB in CB_list:
        f.write(CB + "\n")

#pdb.set_trace()
