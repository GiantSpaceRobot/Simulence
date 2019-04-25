#!/usr/bin/env python

"""
Generate FASTA/FASTQ with realistic coverage distribution
"""

__author__ = "Paul Donovan" 
__maintainer__ = "Paul Donovan"
__email__ = "pauldonovan@rcsi.com"

import sys
import argparse

#Display help and usage
parser = argparse.ArgumentParser(description="Incorrect number of command line arguments")
parser.add_argument('Format (e.g. FASTA/FASTQ)')
parser.add_argument('FASTA input')
parser.add_argument('Read length (e.g "50")')
parser.add_argument('Fold coverage (e.g "100")')
parser.add_argument('FASTA/FASTQ output')
if len(sys.argv[1:]) == 0:
    parser.print_help()
    parser.exit()
args = parser.parse_args()

#Import libraries
from Bio import SeqIO
import sys
from random import *

#Define parameters
readLength = int(sys.argv[3])
foldChange = float(sys.argv[4])
binaryList = [0,1]
nucleotides = "AGTC"

#Define functions
def repeat_to_length(string_to_expand, length):
   return (string_to_expand * ((length/len(string_to_expand))+1))[:length]

def mutator(my_read):
    #Introduce mutations/INDELs into reads (error profiles from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4787001/)
    readString = ""
    for nucleotide in my_read:

        ### Mutation:
        mutateFloat = random()
        if mutateFloat < 0.0035:
            new_nuc = nucleotide
            while str(new_nuc) == str(nucleotide):
                new_nuc = choice(nucleotides)
            nucleotide = new_nuc
        else:
            pass
        
        ### INDELs:
        indelFloat = random()
        if indelFloat < 3.95e-6:   # equivalent to 3.95 x 10^-6
            if (choice(binaryList) == 1): # Flip a coin for heads/tails. If heads, insertion. If tails, deletion.   
                insertion = choice(nucleotides)
                nucleotide = nucleotide + insertion
            else: 
                nucleotide = ""
        
        ### Write (possibly) modified read to string
        readString = readString + nucleotide
    return readString
 
def read_generator(seq, readLen, foldChng):
    #Divide seq length by fold change
    seqLen = len(seq)
    readsFor1X = float(seqLen/readLen) #Number of reads required to get 1X coverage
    readsForFoldChng = int(readsFor1X*foldChng) #Number of reads required to get given fold change coverage
    read_list = list()
    for n in range(int(foldChng)):
        seq1 = seq
        for i in range(int(readsFor1X)):     ### Randomly select sequence from either end of given seq
            if (choice(binaryList) == 1):
                new_read = seq1[0:readLen]
                seq1 = seq1[readLen:]
            else:
                new_read = seq1[-readLen:]
                seq1 = seq1[:-readLen]
            #mutate read
            mutator(new_read)
            read_list.append(new_read)    
    return read_list

fasta_sequences = SeqIO.parse(open(sys.argv[2]),'fasta')
out_file = open(sys.argv[5], "w")
for fasta in fasta_sequences:
    name, sequence, description = fasta.id, str(fasta.seq), str(fasta.description)
    reads = read_generator(sequence, readLength, foldChange)
    count = 1
    for new_seq in reads:
        new_name = name + "_simread" + str(count)
        if sys.argv[1] == "FASTA":
            out_file.write(">" + new_name + "\n" + new_seq + "\n") #FASTA output
        elif sys.argv[1] == "FASTQ":
            quality_score = repeat_to_length("I", readLength)  # Using Q-score of 40 (I) for simulated FASTQ file
            out_file.write("@" + new_name + "\n" + new_seq + "\n+\n" + quality_score + "\n") #FASTQ output
        count = count + 1
out_file.close()

