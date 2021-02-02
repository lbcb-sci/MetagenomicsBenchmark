import sys
import re
import os
from os.path import isfile, join
from Bio import SeqIO

def main_func(dataset, database, genome_sizes_filename, reads_sizes, root_cleaned, root_abundances, names_lines):
    tools = ["truth", "metamaps", "kraken", "clark", "centrifuge", "megan"]

    print("main_func")

    genome_sizes_file = open(genome_sizes_filename, "r")
    genome_sizes_lines = genome_sizes_file.readlines()

    genome_sizes = {}
    for line in genome_sizes_lines:
        parts = re.split(r'\t+', line.strip())
        genome_sizes[parts[0].strip()] = parts[1].strip()


    print("genome_sizes")

    names = {}
    for line in names_lines:
        parts = re.split(r'\|+', line.strip())
        if parts[3].strip() == "scientific name":
            names[parts[0].strip()] = parts[1].strip()

    print("names")

    # print(path_to_dataset)
    # reads_sizes = {}
    # countt = 0
    # for record in SeqIO.parse(path_to_dataset, "fastq"):
    #     if countt % 1000 == 0:
    #         print(countt)
    #     countt += 1
    #     reads_sizes[record.id] = len(record.seq)

    print("reads_sizes")

    results = {}

    print("going tools")

    for tool in tools:
        filename = root_cleaned + "/" + tool + "/" + database + "_" + dataset + ".f2"
        if tool == "truth":
            filename = "truth/" + database + "_" + dataset
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

    print("going results")

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

    filename = root_abundances + "/" + database + "_" + dataset + ".ab"
    file_write = open(filename, "w") 


    print("writing")

    for res in analysis:
        file_write.write(res + "\t" + names[res] + "\t" + str(float(analysis[res][0])) + "\t" + str(float(analysis[res][1])) + "\t" + str(float(analysis[res][2])) + "\t" + str(float(analysis[res][3])) + "\t" + str(float(analysis[res][4])) + "\t" + str(float(analysis[res][5])) + "\n")

    file_write.close()


if __name__ == '__main__':
    if len(sys.argv) < 7:
        print("Script requires 8 arguments")
        exit()
    names_file = open(sys.argv[7], "r")
    names_lines = names_file.readlines()

    read_sizes = {}
    for record in SeqIO.parse(sys.argv[4], sys.argv[8]):
        reads_sizes[record.id] = len(record.seq)
    print("Read reads: " + str(sys.argv[4]))

    main_func(sys.argv[1], sys.argv[2], sys.argv[3], read_sizes, sys.argv[5], sys.argv[6], names_lines)







