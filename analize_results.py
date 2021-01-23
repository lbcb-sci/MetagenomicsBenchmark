import sys
import re
import os

def main_func(root_cleaned, root_reports, dataset, database, names_lines):
    tools = ["truth", "metamaps", "kraken", "clark", "centrifuge", "megan"]

    taxonomy = {}
    for line in names_lines:
        parts = re.split(r'\|+', line.strip())
        if parts[3].strip() == "scientific name":
            taxonomy[parts[0].strip()] = parts[1].strip()

    results = {}

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
                if tax_id in tool_results:
                    tool_results[tax_id] += percentage
                else:
                    tool_results[tax_id] = percentage
        results[tool] = tool_results


    analysis = {}

    for elem in results["truth"]:
        report_list = []
        report_list.append(results["truth"][elem])
        if elem in results["kraken"]:
            report_list.append(results["kraken"][elem])
        else:
            report_list.append(0)
        if elem in results["centrifuge"]:
            report_list.append(results["centrifuge"][elem])
        else:
            report_list.append(0)
        if elem in results["clark"]:
            report_list.append(results["clark"][elem])
        else:
            report_list.append(0)
        if elem in results["metamaps"]:
            report_list.append(results["metamaps"][elem])
        else:
            report_list.append(0)
        if elem in results["megan"]:
            report_list.append(results["megan"][elem])
        else:
            report_list.append(0)
        analysis[elem] = report_list


    for elem in results["kraken"]:
        report_list = []
        if elem in results["truth"]:
            report_list.append(results["truth"][elem])
        else:
            report_list.append(0)
        report_list.append(results["kraken"][elem])
        if elem in results["centrifuge"]:
            report_list.append(results["centrifuge"][elem])
        else:
            report_list.append(0)
        if elem in results["clark"]:
            report_list.append(results["clark"][elem])
        else:
            report_list.append(0)
        if elem in results["metamaps"]:
            report_list.append(results["metamaps"][elem])
        else:
            report_list.append(0)
        if elem in results["megan"]:
            report_list.append(results["megan"][elem])
        else:
            report_list.append(0)
        analysis[elem] = report_list

    for elem in results["centrifuge"]:
        report_list = []
        if elem in results["truth"]:
            report_list.append(results["truth"][elem])
        else:
            report_list.append(0)
        if elem in results["kraken"]:
            report_list.append(results["kraken"][elem])
        else:
            report_list.append(0)
        report_list.append(results["centrifuge"][elem])
        if elem in results["clark"]:
            report_list.append(results["clark"][elem])
        else:
            report_list.append(0)
        if elem in results["metamaps"]:
            report_list.append(results["metamaps"][elem])
        else:
            report_list.append(0)
        if elem in results["megan"]:
            report_list.append(results["megan"][elem])
        else:
            report_list.append(0)
        analysis[elem] = report_list

    for elem in results["clark"]:
        report_list = []
        if elem in results["truth"]:
            report_list.append(results["truth"][elem])
        else:
            report_list.append(0)
        if elem in results["kraken"]:
            report_list.append(results["kraken"][elem])
        else:
            report_list.append(0)
        if elem in results["centrifuge"]:
            report_list.append(results["centrifuge"][elem])
        else:
            report_list.append(0)
        report_list.append(results["clark"][elem])
        if elem in results["metamaps"]:
            report_list.append(results["metamaps"][elem])
        else:
            report_list.append(0)
        if elem in results["megan"]:
            report_list.append(results["megan"][elem])
        else:
            report_list.append(0)
        analysis[elem] = report_list


    for elem in results["metamaps"]:
        report_list = []
        if elem in results["truth"]:
            report_list.append(results["truth"][elem])
        else:
            report_list.append(0)
        if elem in results["kraken"]:
            report_list.append(results["kraken"][elem])
        else:
            report_list.append(0)
        if elem in results["centrifuge"]:
            report_list.append(results["centrifuge"][elem])
        else:
            report_list.append(0)
        if elem in results["clark"]:
            report_list.append(results["clark"][elem])
        else:
            report_list.append(0)
        report_list.append(results["metamaps"][elem])
        if elem in results["megan"]:
            report_list.append(results["megan"][elem])
        else:
            report_list.append(0)
        analysis[elem] = report_list

    for elem in results["megan"]:
        report_list = []
        if elem in results["truth"]:
            report_list.append(results["truth"][elem])
        else:
            report_list.append(0)
        if elem in results["kraken"]:
            report_list.append(results["kraken"][elem])
        else:
            report_list.append(0)
        if elem in results["centrifuge"]:
            report_list.append(results["centrifuge"][elem])
        else:
            report_list.append(0)
        if elem in results["clark"]:
            report_list.append(results["clark"][elem])
        else:
            report_list.append(0)
        if elem in results["metamaps"]:
            report_list.append(results["metamaps"][elem])
        else:
            report_list.append(0)
        report_list.append(results["megan"][elem])
        analysis[elem] = report_list

    filename = root_reports + database + "_" + dataset + ".report"
    file_write = open(filename, "w") 

    for res in analysis:
        file_write.write(res + "\t" + taxonomy[res] + "\t" + str(int(analysis[res][0])) + "\t" + str(int(analysis[res][1])) + "\t" + str(int(analysis[res][2])) + "\t" + str(int(analysis[res][3])) + "\t" + str(int(analysis[res][4])) + "\t" + str(int(analysis[res][5])) + "\n")

    file_write.close()

if __name__ == '__main__':
    if len(sys.argv) < 6:
        print("Script requires 5 arguments")
        exit()
    names_file = open(sys.argv[5], "r")
    names_lines = names_file.readlines()
    main_func(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], names_lines)




