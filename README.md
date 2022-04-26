# eDNA-Species-Identification

This project takes a folder containing barcode sequences, trims each sequence, and combines the sequences under the same barcode

Requirements:

  Pysam
  
  BioPython

Steps:
  1. run main.py
  2. enter the folder containing the barcode## folders

The program will convert from FASTq to FASTA, run it through BLASTn and output the results to XML files
