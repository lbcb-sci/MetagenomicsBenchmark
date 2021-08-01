import sys
import re
import os

def main_func(root_cleaned, root_reports, dataset, database, names_lines, missing_lines, truth_path, target_rank):
    tools = ["truth", "kraken", "centrifuge", "clark", "metamaps", "megan", "minimapA", "minimapM", "ram", "clark-s"]

    missing = {}

    for line in missing_lines:
        parts = re.split(r'\t+', line.strip())
        missing[parts[0].strip()] = parts[1].strip()

    taxonomy = {}
    for line in names_lines:
        parts = re.split(r'\|+', line.strip())
        if parts[3].strip() == "scientific name":
            taxonomy[parts[0].strip()] = parts[1].strip()

    results = {}

    filename = truth_path + "/" + database + "_" + dataset
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
            continue;
            filename = truth_path + "/" + database + "_" + dataset
        file = open(filename, "r") 
        lines = file.readlines()
        
        tool_mappings = {}
        false_positive = 0
        not_in_truth = 0
        total = 0
        for line in lines:
            parts = re.split(r'\t+', line.strip())
            read_id = parts[0].strip()
            tax_id = parts[1].strip()
            percentage = float(parts[2].strip())
            rank = parts[3].strip()
            if rank == target_rank:
                if read_id in truth_res:
                    tool_mappings[read_id] = tax_id
                else:
                    not_in_truth += 1
        for read_id in truth_res:
            tax_id = truth_res[read_id]
            if tax_id in missing:
                total += 1
                if read_id not in tool_mappings:
                    false_positive += 1

        results[tool] = (false_positive, not_in_truth)

    analysis = {}

    filename = root_reports + database + "_" + dataset + ".report_negative"
    file_write = open(filename, "w") 

    for toool in tools:
        if toool == "truth":
            continue;
        (false_positive, not_in_truth, total) = results[toool]
        file_write.write(str(false_positive) + "\t" + str(not_in_truth) + "\t" + str(total) + "\n")

    file_write.close()

if __name__ == '__main__':

    missing_file_path = "data/missing_species.txt"
    # missing_file_path = "data/missing_genus.txt"
    database = "human2"
    root_cleaned = "cleaned_results/"
    # root_cleaned = "cleaned_results_genus/"
    root_reports = "true_negatives/"
    # root_reports = "true_negatives_genus/"
    truth_path = "truth/"
    # truth_path = "truth_genus/"
    target_rank = "species"
    # target_rank = "genus"

    names_file = open("names.dmp", "r")
    names_lines = names_file.readlines()

    missing_file = open(missing_file_path, "r")
    missing_lines = missing_file.readlines()

    for num in range(8, 9):
        dataset = str(num)
        print("Analysing: " + str("human2") + " - " + dataset)
        main_func(root_cleaned, root_reports, dataset, database, names_lines, missing_lines, truth_path, target_rank)




