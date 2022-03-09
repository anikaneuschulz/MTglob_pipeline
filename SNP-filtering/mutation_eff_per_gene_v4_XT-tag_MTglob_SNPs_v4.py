import pysam
import os
import sys
import pdb
import csv
import copy

#this script tells us, what percentage of A C T or G was mutated to what on a per-gene-basis
#it works by calculationg the mutation rate for each read, sums over all reads and divides by the number of reads for each gene
#it does work on a bulk-basis

#it now does work with global MT tags
# it also does not count multimappers anymore, which was a bigger issue then expected!

# it now also can take a blacklist of SNP positions

if len(sys.argv) == 6:
    ifn = sys.argv[1]
    qual_thr = int(sys.argv[2])
    snp_file = sys.argv[3]
    cutoff_rate = int(sys.argv[4])/100
    coverage_cutoff = int(sys.argv[5])
else:
    print("Please provide me with an input .bam file to analyze, a minimum quality score for mutations, as a file that has known SNP positions, the desired cutoff rate in % and the minimum coverage to consider a position valid.")
    sys.exit()

print("This is the SNP aware version.")

scriptpath = os.path.abspath(sys.argv[0])
scriptdir = "/".join(scriptpath.split("/")[:-1])

#chr_locations = scriptdir + "/chrNameLength.txt"
#chr_locations = scriptdir + "/testchromosome.txt"

#chr_locations = open(chr_locations, "r")

## read in SNP file and create a set
snp_set = []
with open(snp_file, "r") as snps:
    next(snps)
    for line in snps:
        lin = line.rstrip().split("\t")
        #allowing only those positions in that have a TC rate over the desired cutoff rate
        if int(lin[1]) < coverage_cutoff:
            snp_set.append(lin[0])
        else:
            if int(lin[2]) > 0:
                if int(lin[2]) / int(lin[1]) >= cutoff_rate:
                    snp_set.append(lin[0])
                    #pdb.set_trace()
    #to speed up membership tests
    snp_set = set(snp_set)

samfile = pysam.AlignmentFile(ifn, "rb")

tab = str.maketrans("ACTG", "TGAC")
mutationsA = ["AC", "AG", "AT"]
mutationsC = ["CA", "CG", "CT"]
mutationsG = ["GA", "GC", "GT"]
mutationsT = ["TA", "TC", "TG"]

mutation_eff_dict = {}
reads_processed = 0
#sys.stdout.write(str(reads_processed) + " reads processed")


#for line in chr_locations:
#    chr_location = line.rstrip().split("\t")

    #this for now ignores mapping points/PCR duplicates and we know it...
for read in samfile.fetch():
    reads_processed+=1
    if reads_processed % 10000 == 0:
        sys.stdout.write("\r" + str(reads_processed) + " reads processed.")
    try:
        # no InDels allowed, also only stuff that maps to a gene in a primary way
        if read.get_tag("XT") and (read.flag == 0 or read.flag == 16) and not any(x in read.cigarstring for x in ["I", "D"]):
            gene = read.get_tag("XT")
            sequence = read.seq
            if read.is_reverse:
                sequence = sequence.translate(tab)[::1]
            cont_A = sequence.count("A")
            cont_C = sequence.count("C")
            cont_G = sequence.count("G")
            cont_T = sequence.count("T")
            mutations = read.get_tag("MT")
            if read.is_reverse:
                mutations = mutations.translate(tab)[::1]
            mutations = mutations.split("_")
            mut_counter = {}
            mut_eff = {}
            mut_counter["AC"] = 0
            mut_counter["AG"] = 0
            mut_counter["AT"] = 0
            mut_counter["CA"] = 0
            mut_counter["CG"] = 0
            mut_counter["CT"] = 0
            mut_counter["GA"] = 0
            mut_counter["GC"] = 0
            mut_counter["GT"] = 0
            mut_counter["TA"] = 0
            mut_counter["TC"] = 0
            mut_counter["TG"] = 0
            for mutation in mutations:
                split_mutation = mutation.split("-")
                if mutation == "0-XX" or mutation == '0-ID':
                    position = 0
                    index = "None"
                else:
                    # InDeld reads are still not allowed to play
                    #if not any(x in read.cigarstring for x in ["I", "D", "N", "S", "H", "P", "X", "B"]) and int(split_mutation[2]) >= qual_thr and not split_mutation[0] in snp_set:
                    if int(split_mutation[2]) >= qual_thr and not split_mutation[0] in snp_set:
                        try:
                            mut_counter[split_mutation[1]]+= 1
                        except KeyError:
                            pass
                            #in case anything containing N reaches high quality, XX or ID

            for mutation in mutationsT:
                try:
                    mut_eff[mutation] = mut_counter[mutation] / cont_T
                except ZeroDivisionError:
                    mut_eff[mutation] = 0
            for mutation in mutationsA:
                try:
                    mut_eff[mutation] = mut_counter[mutation] / cont_A
                except ZeroDivisionError:
                    mut_eff[mutation] = 0
            for mutation in mutationsC:
                try:
                    mut_eff[mutation] = mut_counter[mutation] / cont_C
                except ZeroDivisionError:
                    mut_eff[mutation] = 0
            for mutation in mutationsG:
                try:
                    mut_eff[mutation] = mut_counter[mutation] / cont_G
                except ZeroDivisionError:
                    mut_eff[mutation] = 0


            if gene not in mutation_eff_dict:
                mutation_eff_dict[gene] = []
                mutation_eff_dict[gene].append(copy.deepcopy(mut_eff))
            else:
                mutation_eff_dict[gene].append(copy.deepcopy(mut_eff))
            #if len(mutation_eff_dict) > 5:
            #    pdb.set_trace()

    except KeyError:
        pass
        #stuff that did not map to gene
print("\n")

ofn = ifn[:-4] + "_gene_specific_mutation_rates_Q" + str(qual_thr) + "_SNP-corrected_" + str(cutoff_rate*100) + "perc_cov_" + str(coverage_cutoff) + ".csv"
with open(ofn, "w") as f:
    f.write("gene\tAC\tAG\tAT\tCA\tCG\tCT\tGA\tGC\tGT\tTA\tTC\tTG\treadcount\n")
    for gene in mutation_eff_dict:
        AC = 0
        AG = 0
        AT = 0
        CA = 0
        CG = 0
        CT = 0
        GA = 0
        GC = 0
        GT = 0
        TA = 0
        TC = 0
        TG = 0
        readcount = len(mutation_eff_dict[gene])
        for entry in mutation_eff_dict[gene]:
            AC = AC + entry["AC"]
            AG = AG + entry["AG"]
            AT = AT + entry["AT"]
            CA = CA + entry["CA"]
            CG = CG + entry["CG"]
            CT = CT + entry["CT"]
            GA = GA + entry["GA"]
            GC = GC + entry["GC"]
            GT = GT + entry["GT"]
            TA = TA + entry["TA"]
            TC = TC + entry["TC"]
            TG = TG + entry["TG"]

        AC = AC / len(mutation_eff_dict[gene])
        AG = AG / len(mutation_eff_dict[gene])
        AT = AT / len(mutation_eff_dict[gene])
        CA = CA / len(mutation_eff_dict[gene])
        CG = CG / len(mutation_eff_dict[gene])
        CT = CT / len(mutation_eff_dict[gene])
        GA = GA / len(mutation_eff_dict[gene])
        GC = GC / len(mutation_eff_dict[gene])
        GT = GT / len(mutation_eff_dict[gene])
        TA = TA / len(mutation_eff_dict[gene])
        TC = TC / len(mutation_eff_dict[gene])
        TG = TG / len(mutation_eff_dict[gene])

        f.write(gene + "\t" + str(AC) + "\t" +str(AG) + "\t" +str(AT) + "\t" +str(CA) + "\t" +str(CG) + "\t" +str(CT) + "\t" +str(GA) + "\t" +str(GC) + "\t" +str(GT) + "\t" +str(TA) + "\t" +str(TC) + "\t" + str(TG) + "\t" + str(readcount) +"\n")







#pdb.set_trace()
