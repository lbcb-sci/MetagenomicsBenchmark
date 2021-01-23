import sys
import re
import os
from os import result_listdir
from os.path import isfile, join
from Bio import SeqIO

def main_func(names_lines, dataset, database, genome_sizes_filename, path_to_dataset, tool, root_cleaned):
    tools = ["truth", "metamaps", "kraken", "clark", "centrifuge", "megan"]

    genome_sizes_file = open(genome_sizes_filename, "r")
    genome_sizes_lines = genome_sizes_file.readlines()

    genome_sizes = {}
    for line in genome_sizes_lines:
        parts = re.split(r'\t+', line.strip())
        genome_sizes[parts[0].strip()] = parts[1].strip()

    names = {}
    for line in names_lines:
        parts = re.split(r'\|+', line.strip())
        if parts[3].strip() == "scientific name":
            names[parts[0].strip()] = parts[1].strip()

    reads_sizes = {}
    for record in SeqIO.parse(path_to_dataset, "fasta"):
        reads_sizes[record.id] = len(record.seq)


    results = {}

    for tool in tools:
        filename = root_cleaned + "/" + tool + "/" + database + "_" + dataset + ".f2"
        file_read = open(filename, "r") 

        tool_result = {}

        lines = file_read.readlines()
        for line in lines:
            parts = re.split(r'\t+', line.strip())
            read_id = parts[0].strip()
            tax_id = parts[1].strip()
            percentage = float(parts[2].strip())
            rank = parts[3].strip()
            if rank == "species":
                if read_id in reads_sizes:
                    if tax_id in tool_result:
                        tool_result[tax_id] += (percentage * reads_sizes[read_id])
                    else:
                        tool_result[tax_id] = (percentage * reads_sizes[read_id])

        tool_summary = {}
        for key in tool_result:
            if key in genome_sizes:
                tool_summary[key] = float(tool_result[key]) / float(genome_sizes[key])
        results[tool] = tool_summary

    analysis = {}

    for elem in results["truth"]:
        result_list = []
        result_list.append(results["truth"][elem])
        if elem in results["kraken"]:
            result_list.append(results["kraken"][elem])
        else:
            result_list.append(0)
        if elem in results["centrifuge"]:
            result_list.append(results["centrifuge"][elem])
        else:
            result_list.append(0)
        if elem in results["clark"]:
            result_list.append(results["clark"][elem])
        else:
            result_list.append(0)
        if elem in results["metamaps"]:
            result_list.append(results["metamaps"][elem])
        else:
            result_list.append(0)
        if elem in results["megan"]:
            result_list.append(results["megan"][elem])
        else:
            result_list.append(0)
        analysis[elem] = result_list


    for elem in results["kraken"]:
        result_list = []
        if elem in results["truth"]:
            result_list.append(results["truth"][elem])
        else:
            result_list.append(0)
        result_list.append(results["kraken"][elem])
        if elem in results["centrifuge"]:
            result_list.append(results["centrifuge"][elem])
        else:
            result_list.append(0)
        if elem in results["clark"]:
            result_list.append(results["clark"][elem])
        else:
            result_list.append(0)
        if elem in results["metamaps"]:
            result_list.append(results["metamaps"][elem])
        else:
            result_list.append(0)
        if elem in results["megan"]:
            result_list.append(results["megan"][elem])
        else:
            result_list.append(0)
        analysis[elem] = result_list

    for elem in results["centrifuge"]:
        result_list = []
        if elem in results["truth"]:
            result_list.append(results["truth"][elem])
        else:
            result_list.append(0)
        if elem in results["kraken"]:
            result_list.append(results["kraken"][elem])
        else:
            result_list.append(0)
        result_list.append(results["centrifuge"][elem])
        if elem in results["clark"]:
            result_list.append(results["clark"][elem])
        else:
            result_list.append(0)
        if elem in results["metamaps"]:
            result_list.append(results["metamaps"][elem])
        else:
            result_list.append(0)
        if elem in results["megan"]:
            result_list.append(results["megan"][elem])
        else:
            result_list.append(0)
        analysis[elem] = result_list

    for elem in results["clark"]:
        result_list = []
        if elem in results["truth"]:
            result_list.append(results["truth"][elem])
        else:
            result_list.append(0)
        if elem in results["kraken"]:
            result_list.append(results["kraken"][elem])
        else:
            result_list.append(0)
        if elem in results["centrifuge"]:
            result_list.append(results["centrifuge"][elem])
        else:
            result_list.append(0)
        result_list.append(results["clark"][elem])
        if elem in results["metamaps"]:
            result_list.append(results["metamaps"][elem])
        else:
            result_list.append(0)
        if elem in results["megan"]:
            result_list.append(results["megan"][elem])
        else:
            result_list.append(0)
        analysis[elem] = result_list


    for elem in results["metamaps"]:
        result_list = []
        if elem in results["truth"]:
            result_list.append(results["truth"][elem])
        else:
            result_list.append(0)
        if elem in results["kraken"]:
            result_list.append(results["kraken"][elem])
        else:
            result_list.append(0)
        if elem in results["centrifuge"]:
            result_list.append(results["centrifuge"][elem])
        else:
            result_list.append(0)
        if elem in results["clark"]:
            result_list.append(results["clark"][elem])
        else:
            result_list.append(0)
        result_list.append(results["metamaps"][elem])
        if elem in results["megan"]:
            result_list.append(results["megan"][elem])
        else:
            result_list.append(0)
        analysis[elem] = result_list

    for elem in results["megan"]:
        result_list = []
        if elem in results["truth"]:
            result_list.append(results["truth"][elem])
        else:
            result_list.append(0)
        if elem in results["kraken"]:
            result_list.append(results["kraken"][elem])
        else:
            result_list.append(0)
        if elem in results["centrifuge"]:
            result_list.append(results["centrifuge"][elem])
        else:
            result_list.append(0)
        if elem in results["clark"]:
            result_list.append(results["clark"][elem])
        else:
            result_list.append(0)
        if elem in results["metamaps"]:
            result_list.append(results["metamaps"][elem])
        else:
            result_list.append(0)
        result_list.append(results["megan"][elem])
        analysis[elem] = result_list


    filename = "abundances" + dataset + "_" + database + ".txt"
    file_write = open(filename, "w") 

    for res in analysis:
        file_write.write(res + "\t" + names[res] + "\t" + str(float(analysis[res][0])) + "\t" + str(float(analysis[res][1])) + "\t" + str(float(analysis[res][2])) + "\t" + str(float(analysis[res][3])) + "\t" + str(float(analysis[res][4])) + "\t" + str(float(analysis[res][5])) + "\n")

    file_write.close()


if __name__ == '__main__':
    if len(sys.argv) < 6:
        print("Script requires 5 arguments")
        exit()
    names_file = open(sys.argv[7], "r")
    names_lines = names_file.readlines()
    main_func(names_lines, sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])







