import sys
import re
import os

def main_func(tool, results_path, nodes_lines, fileout, target_rank):
    results_file = open(results_path, "r") 
    results_lines = results_file.readlines()

    outfile = open(fileout,"w") 

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
                if parent_rank == target_rank:
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

    for line in results_lines:
        centrifuge_results = []

        parts = re.split(r'\t+', line.strip())
        is_mapped = parts[0].strip()
        if is_mapped != "C":
            continue

        read_id = parts[1].strip()
        tax_id = parts[2].strip()
        percentage = parts[3].strip()
        if tool != "metamaps":
            percentage = parts[4].strip()

        if tax_id in parents:
            parent = parents[tax_id]
            rank = ranks[tax_id]
            if rank == target_rank:
                outfile.write(read_id.strip() + "\t" + tax_id.strip() + "\t" + percentage.strip() + "\t" + rank.strip() + "\n")
            else:
                resulting_tax_id = find_resulting_tax_id(parent)
                if resulting_tax_id != 0:
                    outfile.write(read_id.strip() + "\t" + resulting_tax_id.strip() + "\t" + percentage.strip() + "\t" + ranks[resulting_tax_id].strip() + "\n")
                else:
                    outfile.write(read_id.strip() + "\t" + tax_id.strip() + "\t" + percentage.strip() + "\t" + rank.strip() + "\n")

    outfile.close()


if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Script requires 5 arguments")
        exit()
    nodes_file = open(sys.argv[3], "r")
    nodes_lines = nodes_file.readlines()
    main_func(sys.argv[1], sys.argv[2], nodes_lines, sys.argv[4], sys.argv[5])

