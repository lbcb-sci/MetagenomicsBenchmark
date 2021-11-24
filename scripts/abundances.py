import analize_abundances
import os
import sys
from Bio import SeqIO

if __name__ == '__main__':
	databases = ["human2", "custom2"]
	root_cleaned = "cleaned_results"
	roots_abundances = "abundances_reports/"
	start = 1

	datasets_files = {
		"1": "01_Mock_100000-bacteria-l1000-q10.fastq",
		"2": "02_silico-30p-human-70p-bac.fastq",
		"3": "03_silico-10-bacteria-100k-reads.fasta",
		"4": "04_silico-3-euka-bac-100k-reads.fastq",
		"5": "05_human-pathogen.fastq",
		"6": "06_50-bac-100k.fastq",
		"8": "08_negative2_10bac_shuffled_human_20k.fasta",
		"9": "09_zymo_pacbio.fastq",
		"10": "10_zymo_ont.fastq",
		"11": "11_SRR11606871_subsampled.fastq"
	}

	extensions = {
		"1": "fastq",
		"2": "fastq",
		"3": "fasta",
		"4": "fastq",
		"5": "fastq",
		"6": "fastq",
		"8": "fasta",
		"9": "fastq",
		"10": "fastq",
		"11": "fastq"
	}

	list_of_read_sizes = {}
	thresholds = {}

	for key in datasets_files:
		value = datasets_files[key]
		reads_sizes = {}
		path_to_dataset = sys.argv[2] + "/" + value
		extension = extensions[key]
		counter = 0
		read_sizes_to_sort = []
		for record in SeqIO.parse(path_to_dataset, extension):
			if counter % 20000 == 0:
				print("Reading: " + str(counter))
			counter += 1
			reads_sizes[record.id] = len(record.seq)
			# read_sizes_to_sort.append(len(record.seq))
		# read_sizes_to_sort.sort()
		threshold = 0
		# threshold = read_sizes_to_sort[int(len(read_sizes_to_sort) * 0.7)]
		thresholds[key] = threshold
		list_of_read_sizes[key] = reads_sizes
		print("Read reads: " + str(path_to_dataset))

	if os.path.isdir(root_cleaned) == False:
		os.mkdir(root_cleaned)

	if os.path.isdir(roots_abundances) == False:
		os.mkdir(roots_abundances)

	names_file = open("names.dmp", "r")
	names_lines = names_file.readlines()
	
	for database in databases:
		for num in range(start, 12):
			if num == 7:
				continue
			dataset = str(num)
			read_sizes = list_of_read_sizes[dataset]
			threshold = thresholds[dataset]
			print("Analysing abundances: " + str(database) + " - " + dataset)
			analize_abundances.main_func(dataset, database, sys.argv[1], read_sizes, threshold, root_cleaned, roots_abundances, names_lines)

