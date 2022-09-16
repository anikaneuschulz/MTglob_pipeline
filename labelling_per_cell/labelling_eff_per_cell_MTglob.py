import pysam
import sys
import pdb

# this calculates labeling efficiency for each cell

ifn = sys.argv[1]
ofn = ifn[:-4] + "_mutation_rate_per_cell.tsv"

qcut = int(sys.argv[2])

bamfile = pysam.AlignmentFile(ifn, "r")

cell_dict = {}
# structure: total #ACGT, #T>C over qcut, #reads

cell_dict2 = {}
# structure: average #ACGT, T>C rate, #reads

#pdb.set_trace()
reads_processed = 0
print("collecting cells...")
sys.stdout.write(str(reads_processed) + " reads processed")
for read in bamfile:
    reads_processed+=1
    sys.stdout.write("\r" + str(reads_processed) + " reads processed")
    if read.has_tag("CB"):
        CB = read.get_tag("CB")
        if not CB in cell_dict:
            cell_dict[CB] = [0,0,0,0,0,0]
        #base counting
        seq = read.seq
        A = seq.count("A")
        C = seq.count("C")
        G = seq.count("G")
        T = seq.count("T")
        cell_dict[CB][0]+=A
        cell_dict[CB][1]+=C
        cell_dict[CB][2]+=G
        cell_dict[CB][3]+=T
        #mutation counting
        mutations = 0
        if read.is_reverse:
            m = "AG"
        else:
            m = "TC"
        MT = read.get_tag("MT").split("_")
        #pdb.set_trace()
        for mutation in MT:
            mut = mutation.split("-")[1:3]
            if mut[0] == m and int(mut[1]) > qcut:
                mutations+=1
        cell_dict[CB][4]+=mutations
        #correct base counting T to C mutations
        cell_dict[CB][1]-=mutations
        cell_dict[CB][3]+=mutations
        #read counting
        cell_dict[CB][5]+=1

bamfile.close()


print("\ncalculating...")
for cell in cell_dict:
    cell_dict2[cell] = [0,0,0,0,0,0]
    try:
        cell_dict2[cell][0] = cell_dict[cell][0]/cell_dict[cell][5]
    except ZeroDivisionError:
        cell_dict2[cell][0] = 0
    try:
        cell_dict2[cell][1] = cell_dict[cell][1]/cell_dict[cell][5]
    except ZeroDivisionError:
        cell_dict2[cell][1] = 0
    try:
        cell_dict2[cell][2] = cell_dict[cell][2]/cell_dict[cell][5]
    except ZeroDivisionError:
        cell_dict2[cell][2] = 0
    try:
        cell_dict2[cell][3] = cell_dict[cell][3]/cell_dict[cell][5]
    except ZeroDivisionError:
        cell_dict2[cell][3] = 0
    try:
        cell_dict2[cell][4] = cell_dict[cell][4]/cell_dict[cell][3]
    except ZeroDivisionError:
        cell_dict2[cell][4] = 0
    cell_dict2[cell][5] = cell_dict[cell][5]

#pdb.set_trace()
print("writing...")
with open(ofn, "w") as out:
    header = "CB\tA\tC\tG\tT\tTtoC\treadcount\n"
    out.write(header)
    for cell in cell_dict2:
        l = cell + "\t" + str(cell_dict2[cell][0]) + "\t" + str(cell_dict2[cell][1]) + "\t" + str(cell_dict2[cell][2]) + "\t" + str(cell_dict2[cell][3]) + "\t" + str(cell_dict2[cell][4]) + "\t" + str(cell_dict2[cell][5]) + "\n"
        out.write(l)
