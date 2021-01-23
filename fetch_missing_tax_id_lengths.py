import sys
import re
import os
from os import listdir
from os.path import isfile, join
from Bio import SeqIO
import urllib.request
import shutil
import gzip

genome_sizes = {}

genome_sizes_file_1 = open(sys.argv[1], "r")
genome_sizes_1 = genome_sizes_file_1.readlines()

for line in genome_sizes_1:
	parts = re.split(r'\t+', line.strip())
	genome_sizes[parts[0].strip()] = parts[0].strip()

tax_ids_file = open(sys.argv[2], "r")
tad_ids_lines = tax_ids_file.readlines()
tax_ids = {}
for line in tad_ids_lines:
	tax_ids[line.strip()] = line.strip()


summary = {}
summary_ftps = {}
summary_file_1 = open(sys.argv[3], "r")
summary_lines_1 = summary_file_1.readlines()

for line in summary_lines_1:
	parts = re.split(r'\t+', line.strip())
	if line[0] == '#':
		continue

	ftpp = ""
	for elem in parts:
		if elem[0:3] == "ftp":
			ftpp = elem
			break

	if parts[6].strip() in summary:
		summary[parts[6].strip()].append(parts[-1].strip())
	else:
		summary[parts[6].strip()] = [parts[-1].strip()]

	if parts[5].strip() in summary:
		summary[parts[5].strip()].append(parts[-1].strip())
	else:
		summary[parts[5].strip()] = [parts[-1].strip()]


	if parts[6].strip() in summary_ftps:
		summary_ftps[parts[6].strip()].append(ftpp)
	else:
		summary_ftps[parts[6].strip()] = [ftpp]

	if parts[5].strip() in summary_ftps:
		summary_ftps[parts[5].strip()].append(ftpp)
	else:
		summary_ftps[parts[5].strip()] = [ftpp]


outfile = sys.argv[4]
file_write_missing = open(outfile + "_missing", "w")
file_write_exceptions = open(outfile + "_exceptions", "w")

not_found = []

for value in tax_ids:
	if value in genome_sizes:
		continue
	else:
		file_write_missing.write(value + "\n")
		not_found.append(value)

for key in not_found:
	if key in summary:
		continue
	else:
		file_write_exceptions.write(key + "\n")

file_write_missing.close()
file_write_exceptions.close()

file_write_new_lengths = open(outfile + "_missing_genome_sizes", "w")

for tax in not_found:
	if tax in summary_ftps:
		length = 0
		count = 0
		summary = summary_ftps[tax]
		for elem in summary:
			parts = re.split(r'/+', elem.strip())
			url = elem + "/" + parts[-1] + "_genomic.fna.gz"
			# filename = wget.download(url)
			with urllib.request.urlopen(url) as response, open("tmp.fna.gz", 'wb') as out_file:
				data = response.read()
				out_file.write(data)
			print(url)
			with gzip.open("tmp.fna.gz", "rt") as handle:
				count += 1
				for record in SeqIO.parse(handle, "fasta"):
					length += len(record.seq)
		file_write_new_lengths.write(tax + "\t" + str(length / count) + "\n")

file_write_new_lengths.close()

