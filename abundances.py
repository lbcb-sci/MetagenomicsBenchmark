import analize_abundances
import os
import sys
from Bio import SeqIO

if __name__ == '__main__':
	tools = ["kraken", "centrifuge", "clark", "metamaps", "megan"]
	databases = ["default", "custom", "human"]
	technologies = ["pb", "ont"]
	root_cleaned = "tax_cleaned_results/"
	roots_abundances = "abundances_reports"

	datasets_files = {
		"ont_1": "Mock_100000-bacteria-l1000-q10.fastq",
		"ont_2": "silico-30p-human-70p-bac.fastq",
		"pb_1": "silico-10-bacteria-100k-reads.fastq",
		"pb_2": "silico-3-euka-bac-100k-reads.fastq",
		"pb_3": "human-pathogen.fastq",
		"pb_4": "50-bac-100k.fastq"
	}

	list_of_read_sizes = {}


	for key in datasets_files:
		value = datasets_files[key]
		reads_sizes = {}
		path_to_dataset = sys.argv[2] + "/" + value
		counter = 0
		for record in SeqIO.parse(path_to_dataset, "fastq"):
			if counter % 20000 == 0:
				print("Reading: " + str(counter))
			counter += 1
			reads_sizes[record.id] = len(record.seq)
		list_of_read_sizes[key] = reads_sizes
		print("Read reads: " + str(path_to_dataset))

	if os.path.isdir(root_cleaned) == False:
		os.mkdir(root_cleaned)

	if os.path.isdir(roots_abundances) == False:
		os.mkdir(roots_abundances)

	names_file = open("names.dmp", "r")
	names_lines = names_file.readlines()
	
	for database in databases:
		for num in range(1, 5):
			for tech in technologies:
				if num > 2 and tech == "ont":
				# if num > 1 or tech == "ont":
					continue
				dataset = tech + "_" + str(num)
				read_sizes = list_of_read_sizes[dataset]
				print("Analysing abundances: " + str(database) + " - " + dataset)
				analize_abundances.main_func(dataset, database, sys.argv[1], read_sizes, root_cleaned, roots_abundances, names_lines)

