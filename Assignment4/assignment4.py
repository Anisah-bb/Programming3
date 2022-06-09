from Bio import SeqIO
import numpy as np
import argparse as ap
import sys


argparser = ap.ArgumentParser(description="Read kmer size")
argparser.add_argument("-kmers", action="store",
                           dest="kmers", required=True, type=int,
                           help="number of kmers used to build this assembly.")
                           
args = argparser.parse_args()
kmer = str(args.kmers)


# define empty list to 
records = []

# read contig file from stdin
output = sys.stdin.readlines()
for line in output:
    if line.startswith('>'):
        start = line.split('N')[0] # splt and strip away the header
        records.append(start) # append > as separator to the record list
    else:
        records.append(line.rstrip('\n')) # append the sequence lines

# make everything a big string and then split it to a list on >
records = ''.join(records)
sequence = records.split('>')
sequence = sequence[1:]

#get lengths of contigs
contigs_len = []
for record in sequence:
    seq = len(record)
    contigs_len.append(seq)

#sort the contigs lens
sorted_contigs_len =  sorted(contigs_len)[::-1]

#calculate N50
def calculate_N50(list_of_lengths):
    '''funtyion to calculate N50 from the list of lengths
    input: list of lengths
    output: N50'''
    a = 0 
    b = np.inf
    for i in range(len(list_of_lengths)):

        if list_of_lengths[i] < b:
            b = list_of_lengths[i] 


        a +=  list_of_lengths[i]

        if a >= sum(list_of_lengths)/2:
            return b




n50  = calculate_N50(sorted_contigs_len)


sys.stdout.write(f"{kmer},{n50}\n") 