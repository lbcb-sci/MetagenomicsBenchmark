import sys
import re
import os

def main_func(root_cleaned, root_reports, dataset, database, names_lines):
    tools = ["truth", "kraken", "centrifuge", "clark", "metamaps", "megan", "minimap2", "minimap", "ram", "clark-s"]

    taxonomy = {}
    for line in names_lines:
        parts = re.split(r'\|+', line.strip())
        if parts[3].strip() == "scientific name":
            taxonomy[parts[0].strip()] = parts[1].strip()

    results = {}

    for tool in tools:
        filename = root_cleaned + tool + "/" + database + "_" + dataset + ".f2"
        if tool == "truth":
            filename = "truth2/" + database + "_" + dataset
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
                if tax_id in tool_results:
                    tool_results[tax_id] += percentage
                else:
                    tool_results[tax_id] = percentage

        results[tool] = tool_results

    if tool == "truth" and (dataset == "9" or dataset == "10" or dataset == "11" or dataset == "12"):
        filename = "truth2/" + database + "_" + dataset + ".f3"
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

    
    filename = root_reports + database + "_" + dataset + ".report"
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




