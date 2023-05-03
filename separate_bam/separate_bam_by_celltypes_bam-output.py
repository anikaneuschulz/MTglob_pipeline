import pysam
import pandas as pd
import sys
import pdb
from subprocess import check_output
import os

ifn = sys.argv[1]
cb_file = sys.argv[2]

#pdb.set_trace()

bamfile = pysam.AlignmentFile(ifn, mode = "rb")

sys.stdout.write("reading cell identities...\n")
cb_ct = pd.read_csv(cb_file, sep = ",", header = None)
cb_ct.columns = ["cb", "celltype"]
cb_dict=dict(zip(cb_ct["cb"],cb_ct["celltype"]))

sys.stdout.write("estimating file size...\n")
linecount = int(str(check_output(["samtools", "view", "-c", ifn]).rstrip()).split("'")[1])
reads_processed = 0
sys.stdout.write("processing...\n")

os.system("mkdir celltype_split")

# create/open handles
handle_dict = {}
for celltype in set(cb_dict.values()):
    print(celltype)
    handle_dict[celltype] = pysam.AlignmentFile("celltype_split/"+celltype+".bam", "wb", template = bamfile)


# iterate
reads_processed = 0
for read in bamfile.fetch():
    reads_processed+=1
    if reads_processed % 10000 == 0:
        sys.stdout.write("\r" + str(round(reads_processed/linecount*100, 2)) + " %")
    #print(read.to_string())
    if read.has_tag("CB"):
        try:
            barcode = read.get_tag("CB")
            #print(barcode)
            celltype = cb_dict[barcode]
            output_file = handle_dict[celltype]
            output_file.write(read)
        #print(cb_dict[barcode])
        except KeyError:
            #in case the file is not 'actual cells' filtered and barcodes don't exist in the matrix
            pass



# closing handles
for celltype in set(cb_dict.values()):
    print(celltype)
    handle_dict[celltype].close()
