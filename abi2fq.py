#!/home/tom/anaconda3/bin/python

############################################################################
# abi2fq.py
# written by Thomas Nelson
# thomas.c.nelson@gmail.com
# This script requires the Biopython module
# Get the Biopython module here here: http://biopython.org/wiki/Download
# This script will convert a .ab1 (abi file) to a fastq file
############################################################################

import sys
from Bio import SeqIO

#usage statement
if len(sys.argv) < 2:
      print("Usage: python abi2fq.py [ abi file name ]")
      sys.exit()

#open file and get data
handle = open(sys.argv[1], "rb")
for record in SeqIO.parse(handle, "abi"):
    name = record.id
    sequence = list(record.seq)
    phred_quality = list(record.letter_annotations.values())[0]
handle.close()

#convert phred scores to phred+33 encoded characters
#phred+33 encoding (0 - 126): list element=phred score, value=phred+33 ascii encoding
phred_33 = ['!','"','#','$','%','&',"'",'(',')','*','+',',','-','.','/','0','1','2','3','4','5','6','7','8','9',':',';','<','=','>','?','@','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','[',"\\",']','^','_','`','a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q','r','s','t','u','v','w','x','y','z','{','|','}','~']
quality = []
for score in phred_quality:
    quality.append(phred_33[score])

#convert lists to strings for easier formatting
sequence_string = (format(''.join(map(str, sequence))))
quality_string = (format(''.join(map(str, quality))))

#write to fastq file
output_file_name = name + ".fastq"
out_handle = open(output_file_name, "w")
out_handle.write("@"+name+"\n")
out_handle.write(sequence_string+"\n")
out_handle.write("+\n")
out_handle.write(quality_string+"\n")
out_handle.close()
