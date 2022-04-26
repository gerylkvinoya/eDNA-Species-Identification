##########################################################
# CS423 Project
# Spring 2022
# 
# blast.py
# 
# This program takes the output from fastq_to_fasta.py
# and runs BLASTn on it, writing the results/data to
# multiple .xml files
#
# Authors:
#   Geryl Vinoya
#   Charlie Benning
#   Gianni Magliana
##########################################################

from Bio.Blast import NCBIWWW
import os
def blastn(directory):
    for file in os.listdir(directory):
        f = os.path.join(directory, file)

        if os.path.isfile(f):
            if f.endswith('.fasta'):

                #create a new filename for the results
                outputFilename = f.replace(".fasta", "_BLASTn_results.xml")

                #open the fasta file
                sequence_data = open(f).read()

                #perform BLASTn on the sequence
                result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)

                #write the results to the file
                writeFile = open(outputFilename, "w")
                writeFile.write(result_handle.read())
                writeFile.close()

    #Source: https://www.delftstack.com/howto/python/remove-substring-from-a-string-python/