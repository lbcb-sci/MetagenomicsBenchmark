import os
import sys
import re


file = open(str(sys.argv[1]), "r")
lines = file.readlines()
outfile = open(str(sys.argv[2]), "w")

nodes_file = open(sys.argv[3], "r")
nodes_lines = nodes_file.readlines()

parents = {}
ranks = {}

for line in nodes_lines:
	parts = re.split(r'\t\|\t+', line.strip())
	parents[parts[0].strip()] = parts[1].strip()
	ranks[parts[0].strip()] = parts[2].strip()

def find_resulting_tax_id(parent):
	not_found_resulting_tax_id = True
	resulting_tax_id = 0
	while not_found_resulting_tax_id:
		if parent in parents:
			parent_rank = ranks[parent]
			if parent_rank == "species":
				not_found_resulting_tax_id = False
				resulting_tax_id = parent
			else:
				prev_parent = parent
				parent = parents[parent]
				if parent == prev_parent:
					not_found_resulting_tax_id = False
		else:
			not_found_resulting_tax_id = False
	return resulting_tax_id

results = {}
max_values = {}

duplicates = 0

for line in lines:
	parts = re.split(r'\t+', line.strip())
	if parts[1].strip() == "0" or parts[1].strip() == "16":
		read_id = parts[0].strip()
		tax_id_extended = parts[2].strip()
	
		parts2 = re.split(r'\|+', tax_id_extended.strip())
		tax_id = parts2[1].strip()
		parent = find_resulting_tax_id(tax_id)

		nm = parts[13].strip()
		parts3 = re.split(r':+', nm.strip())
		value_cig = int(parts3[-1].strip())

		if read_id in results:
			if value_cig > max_values[read_id]:
				results[read_id] = []
				results[read_id].append(parent)
				max_values[read_id] = value_cig
			elif value_cig == max_values[read_id]:
				results[read_id].append(parent)
		else:
			results[read_id] = []
			results[read_id].append(parent)
			max_values[read_id] = value_cig

for read_id in results:
	parents = results[read_id]
	values = {}
	for parent in parents:
		if parent not in values:
			values[parent] = 0
		values[parent] += 1
	max_parent = ""
	max_value = 0
	for parent in values:
		value = values[parent]
		if value > max_value:
			max_value = value
			max_parent = parent
	outfile.write(str(read_id) + "\t" + str(max_parent) + "\t1.0\tspecies\n")

outfile.close()








