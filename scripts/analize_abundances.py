import sys
import re
import os
from os.path import isfile, join
from Bio import SeqIO

def main_func(dataset, database, genome_sizes_filename, reads_sizes, prag, root_cleaned, root_abundances, names_lines):

    tools = ["truth", "kraken", "centrifuge", "clark", "metamaps", "megan", "minimapA", "minimapM", "ram", "clark-s"]

    genome_sizes_file = open(genome_sizes_filename, "r")
    genome_sizes_lines = genome_sizes_file.readlines()

    genome_sizes = {}
    for line in genome_sizes_lines:
        parts = re.split(r'\t+', line.strip())
        genome_sizes[parts[0].strip()] = parts[1].strip()

    print("genome_sizes")

    taxonomy = {}
    for line in names_lines:
        parts = re.split(r'\|+', line.strip())
        if parts[3].strip() == "scientific name":
            taxonomy[parts[0].strip()] = parts[1].strip()

    print("names")

    results = {}

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
                    if int(reads_sizes[read_id]) > prag:
                        if tax_id in tool_result:
                            tool_result[tax_id] += (percentage * reads_sizes[read_id])
                        else:
                            tool_result[tax_id] = (percentage * reads_sizes[read_id])

        tool_summary = {}
        tool_sum = 0.0
        for key in tool_result:
            if key in genome_sizes:
                val_x = float(tool_result[key]) / (float(genome_sizes[key]) * 1000000)
                tool_sum += val_x
            else:
                print("ERROR ! " + str(key) + "\t" + str(tool))
        for key in tool_result:
            if key in genome_sizes:
                val_x = float(tool_result[key]) / (float(genome_sizes[key]) * 1000000)
                tool_summary[key] = (val_x / tool_sum) * 100.0
        results[tool] = tool_summary

        if tool == "truth" and (dataset == "9" or dataset == "10" or dataset == "11"):
            filename = "truth/" + database + "_" + dataset + ".f3"
            file = open(filename, "r") 
            lines = file.readlines()

            tool_results = {}
            for line in lines:
                parts = re.split(r'\t+', line.strip())
                tax_id = parts[0].strip()
                percentage = float(parts[1].strip())
                if tax_id in tool_results:
                    tool_results[tax_id] += percentage
                else:
                    tool_results[tax_id] = percentage

            results[tool] = tool_results

    analysis = {}

    for toool in tools:
        for elem in results[toool]:
            report_list = []
            for toool_inner in tools:
                if toool_inner == toool:
                    report_list.append(results[toool][elem])
                else:
                    if elem in results[toool_inner]:
                        report_list.append(results[toool_inner][elem])
                    else:
                        report_list.append(0)
            analysis[elem] = report_list

    filename = root_abundances + "/" + database + "_" + dataset + ".ab"
    file_write = open(filename, "w") 

    for res in analysis:
        ime = "x"
        if res in taxonomy:
            ime = taxonomy[res]
        rezulting_string = res + "\t" + ime + "\t"

        for smthng in analysis[res]:
            rezulting_string += str(float(smthng)) + "\t"
        rezulting_string += "\n"
        file_write.write(rezulting_string)

    file_write.close()


if __name__ == '__main__':
    if len(sys.argv) < 8:
        print("Script requires 8 arguments")
        exit()
    names_file = open(sys.argv[7], "r")
    names_lines = names_file.readlines()

    read_sizes = {}
    read_sizes_to_sort = []
    for record in SeqIO.parse(sys.argv[4], sys.argv[8]):
        read_sizes[record.id] = len(record.seq)
        read_sizes_to_sort.append(len(record.seq))

    read_sizes_to_sort.sort()
    prag = read_sizes_to_sort[len(read_sizes_to_sort) * 0.7]
    print("Read reads: " + str(sys.argv[4]))

    main_func(sys.argv[1], sys.argv[2], sys.argv[3], read_sizes, prag, sys.argv[5], sys.argv[6], names_lines)







