import pysam
import re
import pdb
import sys
import os
from subprocess import check_output
#import pandas as pd

## input has to be sorted and indexed
# Author: Anika Neuschulz

# this is finally fixing the issue with InDels that could not be processed properly and accurate SNP filtering will be made possible by having global genomic positions in the MTtag

if len(sys.argv) == 2:
    ifn = sys.argv[1]

else:
    print("Please make my life easier by providing an input file (indexed .bam).")
    sys.exit()



if sys.stdout.isatty():
    # You're running in a real terminal
    inp = input("I will place the output at Stdout - are you sure you want me to do that (Y/N)?\n This might flood your terminal - taking measures to catch the output is advised. Press N to learn more.\n Alternatively press F to funnel into a .sam file ending on _MTglob.sam. If your input is > 100 Mio reads, this may be a bad idea.\n")
    if inp == "F" :
        ofn = ifn[:-4] + "_MTglob.sam"
        funneling = False
    elif inp == "Y" :
        ofn = "-"
        funneling = True
    else:
        print("All right, I'll exit now. \n Useful strategies to avoid flooding are the addition of: \n --------- \n | samtools view -Sb | samtools sort -m 2G > outfile.bam \n or \n > outfile.sam \n ---------")
        sys.exit()
else:
    ofn = "-"
    funneling = True
    # You're being piped or redirected




scriptpath = os.path.abspath(sys.argv[0])
scriptdir = "/".join(scriptpath.split("/")[:-1])





samfile = pysam.AlignmentFile(ifn, "rb")

#this creates a .sam file, that then has to be converted to .bam, sorted and indexed


mappedreads = pysam.AlignmentFile(ofn, "w", template=samfile)


if not funneling:
# for the progress counter
    sys.stdout.write("estimating file size...\n")
    linecount = int(str(check_output(["samtools", "view", "-c", ifn]).rstrip()).split("'")[1])
    reads_processed = 0
    sys.stdout.write("processing...\n")


for read in samfile.fetch():
    if not funneling:
        reads_processed+=1
        if reads_processed % 10000 == 0:
            sys.stdout.write("\r" + str(round(reads_processed/linecount*100, 2)) + " %")

    mutations = ["0-XX"]

    # setting the total content tag, based on the reference , to avoid influence of mutations
    total_content = {'A' : 0, 'C' : 0, 'G' : 0, 'T' : 0}
    ref_sequence = read.get_reference_sequence()
    for base in total_content.keys():
        total_content[base] += ref_sequence.count(base)
        TC_tag = ''.join([''.join(key) + str(total_content[key]) + ';' for key in total_content.keys()])[:-1]

    # setting the MT tag if there are mutations present (otherwise it will stay "0-XX")
    if int(read.get_tag("NM")) > 0:
        read_sequence = read.query_alignment_sequence
        global_positions = read.get_reference_positions()
        contig = read.reference_name
        md = re.findall('(\d+)(\^[A-Za-z]+|[A-Za-z])', read.get_tag("MD"))
        cumsum_start = -1
        mutations = []
            #first MD number does not need to be altered, all others need to have 1 added for position
            #of mutated base, therefore this cheap -1 trick
        for element in md:
            try:
                base = element[1]
                #filtering out deletions marked by ^ in the MD tag
                if not base[0] == "^":
                    #pdb.set_trace()
                    position = int(element[0]) + cumsum_start + 1
                    sampleBase = read_sequence[position]
                    seq_quality = read.query_alignment_qualities[position]
                    # here we are getting the genomic position from the mapping positions, this should be InDel-safe
                    global_position = contig + ":" + str(global_positions[position])
                    mutations.append(str(global_position) + "-" + base + sampleBase + "-" + str(seq_quality))
                    cumsum_start = position
            except IndexError:
                pdb.set_trace()
    #if "^" in mutations:
#        pdb.set_trace()
    # in case the only modification was a deletion
    if mutations == []:
        mutations = ["0-XX"]
    read.tags = read.tags + [('MT','_'.join(mutations))]
    read.set_tag('TC',TC_tag,'Z')
    mappedreads.write(read)

mappedreads.close()
