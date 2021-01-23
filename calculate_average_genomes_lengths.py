import sys
import re
import os
from Bio import SeqIO

database_path = sys.argv[1]
names_path = sys.argv[2]
names_file = open(names_path, "r") 
names_lines = names_file.readlines()
names = {}
for line in names_lines:
    parts = re.split(r'\t\|\t+', line.strip())
    names[parts[1].strip()] = parts[0].strip()

filename_out = sys.argv[3]
file_write = open(filename_out,"w") 

legend = {}
content = []

onlyfiles = [f for f in os.listdir(database_path) if os.path.isfile(os.path.join(database_path, f))]

for filename in onlyfiles:
    fullpath = database_path + "/" + filename
    fullfilename, file_extension = os.path.splitext(fullpath)

    if(file_extension == ".fna"):
        with open(fullpath, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                parts = re.split(r' +', record.description)
                for i in range(1, (len(parts)-2)):
                    for j in range(i+2, len(parts)):
                        subarr = ""
                        for k in range(i, j):
                            subarr = subarr + " " + parts[k]
                        if subarr[0] == ',':
                            subarr = subarr[1:]
                        if subarr[-1] == ',':
                            subarr = subarr[:-1]
                        if subarr.strip() in names:
                            found_tax_id = names[subarr.strip()]
                            if found_tax_id in legend:
                                legend[found_tax_id].append(fullpath)
                            else:
                                legend[found_tax_id] = [fullpath]
                break

for tax_id in legend:
    total_length = 0
    sequence_number = 0
    genomes_list = legend[tax_id]
    for genome_path in genomes_list:
        fullfilename, file_extension = os.path.splitext(genome_path)
        if(file_extension == ".fna"):
            with open(genome_path, "rU") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    total_length += len(record.seq)
                    sequence_number += 1
    average_length = total_length / sequence_number
    file_write.write(tax_id + "\t" + str(average_length) + "\n")
file_write.close()
