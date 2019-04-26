#!/usr/bin/env python

#Copyright 2019 Paul Donovan
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


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
    read_list = list()
    for n in range(int(foldChng)):
        seq1 = seq
        seqLen = len(seq1)
        readsFor1X = float(seqLen)/float(readLen) #Number of reads required to get 1X coverage
        if readsFor1X < 5:
            newReadLen = seqLen/int(readsFor1X)
            #otherReadLen = seqLen - newReadLen
            print (seqLen, readLen, newReadLen)
        else:
            newReadLen = readLen
        readsForFoldChng = int(readsFor1X*foldChng) #Number of reads required to get given fold change coverage
        for i in range(int(readsFor1X)):     ### Randomly select sequence from either end of given seq
            if (choice(binaryList) == 1):
                new_read = seq1[0:newReadLen]
                seq1 = seq1[newReadLen:]
            else:
                new_read = seq1[-newReadLen:]
                seq1 = seq1[:-newReadLen]
            #mutate read
            mutator(new_read)
            read_list.append(new_read)    
    return read_list

fasta_sequences = SeqIO.parse(open(sys.argv[2]),'fasta')
out_file = open(sys.argv[5], "w")
for fasta in fasta_sequences:
    name, sequence, description = fasta.id, str(fasta.seq), str(fasta.description)
    #if len(sequence) < readLength: # If the provided readLength is longer than the current sequence, change provided read Length
    #    readLength = len(sequence)
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

