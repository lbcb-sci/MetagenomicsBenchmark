import sys
import re
import os

def main_func(root_cleaned, root_reports, dataset, database, names_lines):
    tools = ["truth", "kraken", "centrifuge", "clark", "metamaps", "megan", "minimapA", "minimapM", "ram", "clark-s"]

    taxonomy = {}
    for line in names_lines:
        parts = re.split(r'\|+', line.strip())
        if parts[3].strip() == "scientific name":
            taxonomy[parts[0].strip()] = parts[1].strip()

    results = {}

    filename = "truth/" + database + "_" + dataset
    file_read = open(filename, "r") 
    truth_res = {}
    lines = file_read.readlines()
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        read_id = parts[0].strip()
        tax_id = parts[1].strip()
        percentage = float(parts[2].strip())
        rank = parts[3].strip()
        truth_res[read_id] = tax_id

    for tool in tools:
        filename = root_cleaned + tool + "/" + database + "_" + dataset + ".f2"
        if tool == "truth":
            filename = "truth/" + database + "_" + dataset
        file = open(filename, "r") 
        lines = file.readlines()
        
        tool_results = {}
        for line in lines:
            parts = re.split(r'\t+', line.strip())
            read_id = parts[0].strip()
            tax_id = parts[1].strip()
            percentage = float(parts[2].strip())
            rank = parts[3].strip()
            if rank == "species":
                if read_id in truth_res:
                    if truth_res[read_id] == tax_id:
                        if tax_id in tool_results:
                            tool_results[tax_id] += percentage
                        else:
                            tool_results[tax_id] = percentage
        results[tool] = tool_results

    analysis = {}

    for tool in tools:
        for elem in results[tool]:
            report_list = []
            for tool_inner in tools:
                if tool_inner == tool:
                    report_list.append(results[tool][elem])
                else:
                    if elem in results[tool_inner]:
                        report_list.append(results[tool_inner][elem])
                    else:
                        report_list.append(0)
            analysis[elem] = report_list

    filename = root_reports + database + "_" + dataset + ".report_truth"
    file_write = open(filename, "w") 

    for res in analysis:
        ime = "x"
        if res in taxonomy:
            ime = taxonomy[res]
        rezulting_string = res + "\t" + ime + "\t"

        for smthng in analysis[res]:
            rezulting_string += str(int(smthng)) + "\t"
        rezulting_string += "\n"
        file_write.write(rezulting_string)

    file_write.close()

if __name__ == '__main__':
    if len(sys.argv) < 6:
        print("Script requires 5 arguments")
        exit()
    names_file = open(sys.argv[5], "r")
    names_lines = names_file.readlines()
    main_func(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], names_lines)



