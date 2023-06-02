import pysam
import re
import pdb
import string
import sys
import os
import random



if len(sys.argv) == 4:
    ifn = sys.argv[1]
    qcut= int(sys.argv[2])
    MT_separator = sys.argv[3]
else:
    print("Please give me an input file to excercise my reverting abilities on, a minimum quality as well as the MT-tag separator (usually _ , except for the refSeq capable scripts which use §).")
    sys.exit()

ofn = ifn.replace("MTglob", "").replace("MD", "")[:-4] + "ARTIFICIAL_CONTROL.bam"
bamfile = pysam.AlignmentFile(ifn, "rb")
outfile = pysam.AlignmentFile(ofn, "wb", template=bamfile)


## this is where we set up TG and TA counting - the background TC rate is usually right between the two, so we can exploit that to estimate the nessesary mutation rate

mutation_dict = {}
mutation_dict["TA"] = 0
mutation_dict["TG"] = 0

tab = str.maketrans("ACTG", "TGAC")

countA = 0
countT = 0
countC = 0
countG = 0

sys.stdout.write("estimating background mutation rate...\n")
lines_processed=0
for read in bamfile.fetch():
    lines_processed+=1
    if lines_processed % 500000 == 0:
        sys.stdout.write('\r' + str(lines_processed) + ' reads processed.')

    if read.flag == 0 or read.flag == 16:
        sequence = read.seq
        if not read.is_reverse:
            countA+= sequence.count("A")
            countT+= sequence.count("T")
            countC+= sequence.count("C")
            countG+= sequence.count("G")

        else:
            countT+= sequence.count("A")
            countA+= sequence.count("T")
            countG+= sequence.count("C")
            countC+= sequence.count("C")
        mutations = read.get_tag("MT").split("§")
        ## this is where we catch unmutated reads or insertions/deletions
        for mutation in mutations:
            if mutation == "0-XX" or mutation == "0-ID":
                break
            split_mutation = mutation.split("-")
            # here we adujst for the absolute MT tags
            #mapping_start = read.reference_start
            #relative_position = int(split_mutation[0].split(':')[1]) - mapping_start
            try:
                if int(split_mutation[2]) >= qcut:
                    try:
                        if read.is_reverse:
                            split_mutation[1] = split_mutation[1].translate(tab)[::1]
                        mutation_dict[split_mutation[1]]+=1
                    except KeyError:
                        pass # only TA and TG mutations are bring written to dict, hope that saves some time
            except IndexError:
                pdb.set_trace()

bg_mut_rate = (mutation_dict["TA"]/countT + mutation_dict["TG"]/countT) / 2

sys.stdout.write("\n estimated background mutation rate: " + str(round(bg_mut_rate*100, 6)) + "%")

## reverting mutations, keeping the estimated background mutation rate in

sys.stdout.write("\nreverting mutations...\n")

random.seed(0)
###scaling up the bg mutation rate, as the random numger generator only works with integers, also rounding it for closest results
bg_mut_rate_scaled = round(bg_mut_rate*10000, 0)

lines_processed=0
for read in bamfile.fetch():
    lines_processed+=1
    if lines_processed % 500000 == 0:
        sys.stdout.write('\r' + str(lines_processed) + ' reads processed.')
    sequence = read.seq
    mutations = read.get_tag("MT").split("§")
    mapping_start = read.get_reference_positions()[0]
    ## set this switch to see if the sequence was altered and we need to re-write it at the end of the loop
    tampered_with = False
    for mutation in mutations:
        if mutation == "0-XX" or mutation == "0-ID":
            break
        split_mutation = mutation.split("-")
        if read.is_reverse:
            label = "AG"
            reversal = "A"
        else:
            label = "TC"
            reversal = "T"
        if split_mutation[1] == label:
            TC_in_contig_position = int(split_mutation[0].split(":")[1])
            ## shall we revert this one?
            if random.randrange(0, 10000, 1) > bg_mut_rate_scaled :
                # saving the quality string, in case the bug where it disappeared if the sequence was tampered with, still exists
                quality_stashed = read.qual
                # finding where the TC is located in the sequence, massively inefficient. Maybe there is another robust way, but none has shown up so far
                # we are turning the mapping positions in a dictionary and swapping keys and items afterwards to have the genomic positions searchable
                # there HAS to be another way...
                mapping_positions = dict((v,k) for k,v in dict(read.get_aligned_pairs()).items())
                TC_in_sequence_position = mapping_positions[TC_in_contig_position]
                ## a little detour since strings can't be changed in python...
                sequence_1 = list(sequence)
                sequence_1[TC_in_sequence_position] = reversal
                sequence = "".join(sequence_1)
                tampered_with = True
    if tampered_with == True :
        read.seq = sequence
        read.qual = quality_stashed
    ## removing MD, TC and MT tags, as they may not match the updated sequences any more and will need to be re-generated
    ## the CIGAR string should not care (at least the ones generated by STAR don't differentiate between matches and mismatches)
    read.set_tag("MD", None )
    read.set_tag("MT", None )
    read.set_tag("TC", None )
    read.set_tag("NM", None )
    outfile.write(read)

outfile.close()
