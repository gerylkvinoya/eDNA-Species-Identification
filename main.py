##########################################################
# CS423 Project
# Spring 2022
# 
# main.py
# 
# This is the main program to run
#
# Authors:
#   Geryl Vinoya
#   Charlie Benning
#   Gianni Magliana
##########################################################

from pysam import FastxFile
from Bio.Blast import NCBIWWW
import textwrap
import os
import fastq_to_fasta
import blast

def main():
    directory = input("Enter folder containing barcode sequences:\n")

    fastq_to_fasta.searchDir(directory)

    print("FASTq to FASTA Complete. Running BLASTn")

    blast.blastn(directory)

    print("BLASTn complete")

main()

