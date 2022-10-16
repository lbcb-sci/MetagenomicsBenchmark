import sys
import re
import os
from Bio import SeqIO

def find_resulting_tax_id(current_tax, target_rank, taxonomy_tree, ranks):
    not_found_resulting_tax_id = True
    resulting_tax_id = 0
    tax = current_tax
    while not_found_resulting_tax_id:
        if tax in taxonomy_tree:
            parent_rank = ranks[tax]
            if parent_rank == target_rank:
                not_found_resulting_tax_id = False
                resulting_tax_id = tax
            else:
                prev_tax = tax
                tax = taxonomy_tree[tax]
                if tax == prev_tax:
                    not_found_resulting_tax_id = False
        else:
            not_found_resulting_tax_id = False
    return resulting_tax_id

def transform_results_to_target_rank(parsed_rows, taxonomy_tree, ranks, target_rank):
    transformed_rows = {}
    for row in parsed_rows:
        read_id = row.read_id
        tax_id = row.tax_id

        if tax_id in taxonomy_tree:
            rank = ranks[tax_id]
            if rank == target_rank:
                transformed_rows[read_id.strip()] = (tax_id.strip(), rank.strip())
            else:
                parent_tax_id = taxonomy_tree[tax_id]
                resulting_tax_id = find_resulting_tax_id(parent_tax_id, target_rank, taxonomy_tree, ranks)
                if resulting_tax_id != 0:
                    transformed_rows[read_id.strip()] = (resulting_tax_id.strip(), ranks[resulting_tax_id].strip())
                else:
                    transformed_rows[read_id.strip()] = (tax_id.strip(), rank.strip())
    return transformed_rows


class Row:
    def __init__(self, read_id, tax_id):
        self.read_id = read_id
        self.tax_id = tax_id

    def to_string(self):
        return str(self.read_id) + "\t" + str(self.tax_id) + "\n"

def analyse_metamaps(lines):
    read_first_line = False
    rows = []
    for line in lines:
        if read_first_line == False:
            read_first_line = True
            continue;
        parts = re.split(r'\t+', line.strip())
        row = Row(parts[0].strip(), parts[1].strip())
        rows.append(row)
    return rows

def analyse_centrifuge(lines):
    read_first_line = False
    rows = []
    for line in lines:
        if read_first_line == False:
            read_first_line = True
            continue;
        parts = re.split(r'\t+', line.strip())
        row = Row(parts[0].strip(), parts[2].strip())
        rows.append(row)
    return rows

def analyse_clark(lines):
    read_first_line = False
    rows = []
    for line in lines:
        if read_first_line == False:
            read_first_line = True
            continue;
        parts = re.split(r'\,+', line.strip())
        tax_id = parts[2].strip()
        if tax_id == "NA":
            continue
        row = Row(parts[0].strip(), tax_id)
        rows.append(row)
    return rows

def analyse_clark_s(lines):
    read_first_line = False
    rows = []
    for line in lines:
        if read_first_line == False:
            read_first_line = True
            continue;
        parts = re.split(r'\,+', line.strip())
        tax_id = parts[3].strip()
        if tax_id == "NA":
            continue
        row = Row(parts[0].strip(), tax_id)
        rows.append(row)
    return rows

def analyse_kraken(lines):
    rows = []
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        row = Row(parts[1].strip(), parts[2].strip())
        rows.append(row)
    return rows

def analyse_megan(lines):
    rows = []
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        row = Row(parts[0].strip(), parts[1].strip())
        rows.append(row)
    return rows

def analyse_kaiju(lines):
    rows = []
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        if parts[0].strip() == "C":
            row = Row(parts[1].strip(), parts[2].strip())
            rows.append(row)
    return rows

def analyse_truth_kaiju_exception(lines, target_rank, taxonomy_tree, ranks):
    kaiju_results = {}
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        read_id = parts[0].strip()
        tax_id = parts[1].strip()
        rank = parts[3].strip()
        resulting_tax_id = find_resulting_tax_id(tax_id, target_rank, taxonomy_tree, ranks)
        parts_read_id = re.split(r'/+', read_id)
        kaiju_read_id = "_".join(parts_read_id)
        kaiju_results[kaiju_read_id] = (resulting_tax_id, ranks[resulting_tax_id])
    return kaiju_results

def analyse_truth(lines, target_rank, taxonomy_tree, ranks, is_percentage):
    results = {}
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        if is_percentage:
            tax_id = parts[0].strip()
            percentage = parts[1].strip()
            results[tax_id] = percentage
        else:
            read_id = parts[0].strip()
            tax_id = parts[1].strip()
            rank = parts[3].strip()
            resulting_tax_id = find_resulting_tax_id(tax_id, target_rank, taxonomy_tree, ranks)
            results[read_id] = (resulting_tax_id, ranks[resulting_tax_id])
    return results

def analyse_bracken(lines, target_rank, taxonomy_tree, ranks, is_percentage):
    transformed_rows = {}
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        if parts[0].strip() == "name":
            continue
        read_count = float(parts[5].strip())
        abud = float(parts[6].strip())
        tax_id = parts[1].strip()
        if parts[2].strip() == "S" and target_rank == "species": 
            transformed_rows[tax_id] = (str(abud), str(read_count))
        else:
            resulting_tax_id = find_resulting_tax_id(tax_id, target_rank, taxonomy_tree, ranks)
            if resulting_tax_id != 0:
                transformed_rows[resulting_tax_id] = (str(abud), str(read_count))
    return transformed_rows

def analyse_sam_paf(lines, target_rank, taxonomy_tree, ranks, is_sam):
    results = {}
    max_values = {}
    transformed_rows = {}
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        if is_sam:
            if parts[1].strip() != "0" and parts[1].strip() != "16":
                continue
        read_id = parts[0].strip()
        tax_id_extended = parts[2].strip()
        if is_sam == False:
            tax_id_extended = parts[5].strip()
    
        parts2 = re.split(r'\|+', tax_id_extended.strip())
        tax_id = parts2[1].strip()
        resulting_tax_id = find_resulting_tax_id(tax_id, target_rank, taxonomy_tree, ranks)
        if resulting_tax_id == 0:
            resulting_tax_id = tax_id

        value_cig = 0.0
        if is_sam:
            if len(parts) >= 14 :
                nm = parts[13].strip()
                parts3 = re.split(r':+', nm.strip())
                value_cig = int(parts3[-1].strip())
        else:
            length_q = int(parts[3].strip()) - int(parts[2].strip())
            length_t = int(parts[9].strip()) - int(parts[8].strip())
            length = max(length_t, length_q)
            nm = int(parts[9].strip())
            value_cig = ((2.0 * float(length)) * float(nm)) / (float(length) + float(nm)) 

        if read_id in results:
            if value_cig > max_values[read_id]:
                results[read_id] = []
                results[read_id].append(resulting_tax_id)
                max_values[read_id] = value_cig
            elif value_cig == max_values[read_id]:
                results[read_id].append(resulting_tax_id)
        else:
            results[read_id] = []
            results[read_id].append(resulting_tax_id)
            max_values[read_id] = value_cig

    for read_id in results:
        tax_ids = results[read_id]
        values = {}
        for tax_id in tax_ids:
            if tax_id not in values:
                values[tax_id] = 0
            values[tax_id] += 1
        max_tax_id = ""
        max_value = 0
        for tax_id in values:
            value = values[tax_id]
            if value > max_value:
                max_value = value
                max_tax_id = tax_id
        transformed_rows[read_id.strip()] = (max_tax_id, ranks[max_tax_id])
    return transformed_rows


def read_read_sizes(dataset_path, datasets_path):
    reads_sizes = {}
    reads_sizes_kaiju = {}
    path_to_dataset = datasets_path + "/" + dataset_path
    parts_path = re.split(r'\.+', dataset_path.strip())
    extension = parts_path[-1]
    counter = 0
    read_sizes_to_sort = []
    for record in SeqIO.parse(path_to_dataset, extension):
        counter += 1

        parts_read_id = re.split(r'/+', str(record.id).strip())
        read_id_kaiju = "_".join(parts_read_id)
        reads_sizes_kaiju[read_id_kaiju] = len(record.seq)

        reads_sizes[record.id] = len(record.seq)

        read_sizes_to_sort.append(len(record.seq))
    read_sizes_to_sort.sort()
    threshold = read_sizes_to_sort[int(len(read_sizes_to_sort) * 0.7)]
    return (reads_sizes, reads_sizes_kaiju, threshold)

def read_cleaned_results(number_of_dataset, database, target_rank, tool, root_cleaned_results, is_percentage_tool):
    cleaned_filename = root_cleaned_results + "/" + str(tool) + "_" + str(database) + "_" + str(number_of_dataset) + "_" + str(target_rank) + ".f2"
    cleaned_outfile = open(cleaned_filename, "r")
    lines = cleaned_outfile.readlines()
    transformed_rows = {}
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        if is_percentage_tool:
            transformed_rows[parts[0].strip()] = parts[1].strip()
        else:
            transformed_rows[parts[0].strip()] = (parts[1].strip(), parts[2].strip())
    return transformed_rows

def get_truth_tax_id(all_transformed_rows, tool, seq_id, is_kaiju_exception):
    if tool == "kaiju" and is_kaiju_exception:
        parts_read_id = re.split(r'/+', seq_id)
        kaiju_seq_id = "_".join(parts_read_id)
        if kaiju_seq_id in all_transformed_rows["truthK"]:
            return all_transformed_rows["truthK"][seq_id][0]
    else:
        if seq_id in all_transformed_rows["truth"]:
            return all_transformed_rows["truth"][seq_id][0]
    return ""


def generate_read_count_reports(tools, all_transformed_rows, root_reports, database, dataset, target_rank, is_percentage, percentage_tools, is_kaiju_exception, number_of_nodes_files, taxonomy_names_lists, missing_tax_ids_total, is_negative):
    total_tool_results = {}
    total_tax_ids = {}
    total_tool_TP_results = {}
    total_TP_tax_ids = {}

    report_tools = {}
    report_tools_TP = {}

    true_negatives = {}

    for tool in tools:
        if tool == "truth" and is_percentage == True:
            report_tools[tool] = 1
            tool_results = {}
            for tax_id in all_transformed_rows[tool]:
                percentage = float(all_transformed_rows[tool][tax_id])
                if tax_id in tool_results:
                    tool_results[tax_id] += percentage
                else:
                    tool_results[tax_id] = percentage
                total_tax_ids[tax_id] = 1
            total_tool_results[tool] = tool_results
            continue
        if tool in percentage_tools:
            report_tools[tool] = 1
            tool_results = {}
            for this_tax_id in all_transformed_rows[tool]:
                (abundance, read_count) = all_transformed_rows[tool][this_tax_id]
                if this_tax_id not in tool_results:
                    tool_results[this_tax_id] = 0
                total_tax_ids[this_tax_id] = 1
                tool_results[this_tax_id] += float(read_count)
            
            total_tool_results[tool] = tool_results
            continue

        truth_tool_key = "truthK" if (tool == "kaiju" and is_kaiju_exception) else "truth"
        report_tools[tool] = 1
        report_tools_TP[tool] = 1
        tool_results = {}
        tool_TP_results = {}
        tool_mappings = {}
        not_in_truth = 0
        true_negative = 0
        for this_seq_id in all_transformed_rows[tool]:
            (this_tax_id, this_rank) = all_transformed_rows[tool][this_seq_id]
            if this_rank != target_rank:
                continue

            if this_tax_id not in tool_results:
                tool_results[this_tax_id] = 0
            total_tax_ids[this_tax_id] = 1
            tool_results[this_tax_id] += 1

            if this_seq_id in all_transformed_rows[truth_tool_key]:
                tool_mappings[this_seq_id] = 1
            else:
                not_in_truth += 1

            truth_tax_id = get_truth_tax_id(all_transformed_rows, tool, this_seq_id, is_kaiju_exception)
            if truth_tax_id == "":
                continue

            if truth_tax_id == this_tax_id:
                if this_tax_id not in tool_TP_results:
                    tool_TP_results[this_tax_id] = 0
                total_TP_tax_ids[this_tax_id] = 1
                tool_TP_results[this_tax_id] += 1
        
        total_tool_results[tool] = tool_results
        total_tool_TP_results[tool] = tool_TP_results

        for this_seq_id in all_transformed_rows[truth_tool_key]:
            truth_tax_id = all_transformed_rows[truth_tool_key][this_seq_id][0]
            if truth_tax_id in missing_tax_ids_total[tools[tool]][target_rank]:
                if this_seq_id not in tool_mappings:
                    true_negative += 1
        true_negatives[tool] = (true_negative, 20000 - not_in_truth)

    report = generate_report(total_tax_ids, total_tool_results)
    report_TP = generate_report(total_TP_tax_ids, total_tool_TP_results)

    read_count_filename = root_reports + "/" + database + "_" + dataset + "_" + target_rank + "_read_count.csv"
    write_report(report, number_of_nodes_files, taxonomy_names_lists, read_count_filename, report_tools)
    read_count_TP_filename = root_reports + "/" + database + "_" + dataset + "_" + target_rank + "_read_count_TP.csv"
    write_report(report_TP, number_of_nodes_files, taxonomy_names_lists, read_count_TP_filename, report_tools_TP)
    read_count_TN_filename = root_reports + "/" + database + "_" + dataset + "_" + target_rank + "_read_count_TN.csv"
    write_report_TN(true_negatives, number_of_nodes_files, taxonomy_names_lists, read_count_TN_filename, report_tools_TP, is_negative)

def generate_report(total_tax_ids, total_tool_results):
    report = {}
    for this_tax_id in total_tax_ids:
        this_tax_id_results = []
        for tool in total_tool_results:
            if this_tax_id in total_tool_results[tool]:
                this_tax_id_results.append(total_tool_results[tool][this_tax_id])
            else:
                this_tax_id_results.append(0)
        report[this_tax_id] = this_tax_id_results
    return report

def write_report_TN(true_negatives, number_of_nodes_files, taxonomy_names_lists, read_count_filename, report_tools, is_negative):
    read_count_file = open(read_count_filename, "w")
    for tool in report_tools:
        if tool in true_negatives:
            value = true_negatives[tool][1] if is_negative else true_negatives[tool][0]
            read_count_file.write(str(tool) + "\t" + str(value) + "\n")
    read_count_file.close()

def write_report(report, number_of_nodes_files, taxonomy_names_lists, read_count_filename, report_tools):
    read_count_file = open(read_count_filename, "w")
    read_count_file.write("tax id\tname\t")
    for tool in report_tools:
        read_count_file.write(str(tool) + "\t")
    read_count_file.write("\n")
    for this_tax_id in report:
        name = "x"
        rezulting_string = ""
        for names_index in range(0, int(number_of_nodes_files)):
            if this_tax_id in taxonomy_names_lists[names_index]:
                name = taxonomy_names_lists[names_index][this_tax_id]
                break
        rezulting_string = this_tax_id + "\t" + name + "\t"

        for result_value in report[this_tax_id]:
            rezulting_string += str(float(result_value)) + "\t"
        rezulting_string += "\n"
        read_count_file.write(rezulting_string)
    read_count_file.close()

def generate_abundance_reports(tools, all_transformed_rows, root_reports, database, dataset, is_percentage, percentage_tools, is_kaiju_exception, number_of_nodes_files, taxonomy_names_lists, reads_sizes, reads_sizes_kaiju, threshold, genome_sizes):
    results_abundance = {}
    results_read_count = {}
    results_threshold = {}

    total_tax_ids_abundance = {}
    total_tax_ids_threshold = {}

    report_tools_abundance = {}
    report_tools_read_count = {}

    for tool in tools:
        target_reads_sizes = reads_sizes
        if tool == "kaiju":
            target_reads_sizes = reads_sizes_kaiju
        if tool in percentage_tools:
            report_tools_abundance[tool] = 1
            tool_results = {}
            tool_summary = {}
            tool_sum = 0.0
            tool_results_read_count = {}
            tool_summary_read_count = {}
            tool_sum_read_count = 0.0
            report_tools_read_count[tool] = 1
            for tax_id in all_transformed_rows[tool]:
                (percentage, read_cnt_perc) = all_transformed_rows[tool][tax_id]
                if tax_id in genome_sizes:
                    abundance = float(percentage) / (float(genome_sizes[tax_id]) * 1000000)
                    tool_results[tax_id] = abundance
                    tool_sum += abundance
                    tool_results_read_count[tax_id] = float(percentage)
                    tool_sum_read_count += float(percentage)
            for tax_id in tool_results:
                percentage = tool_results[tax_id]
                tool_summary[tax_id] = (percentage / tool_sum) * 100
                total_tax_ids_abundance[tax_id] = 1
            results_abundance[tool] = tool_summary
            for tax_id in tool_results_read_count:
                percentage = tool_results_read_count[tax_id]
                tool_summary_read_count[tax_id] = (percentage / tool_sum_read_count) * 100
                total_tax_ids_abundance[tax_id] = 1
            results_read_count[tool] = tool_summary_read_count
        elif tool == "truth" and is_percentage: 
            report_tools_abundance[tool] = 1
            report_tools_read_count[tool] = 1  
            tool_results = {}
            for tax_id in all_transformed_rows[tool]:
                percentage = float(all_transformed_rows[tool][tax_id])
                if tax_id in tool_results:
                    tool_results[tax_id] += percentage
                else:
                    tool_results[tax_id] = percentage
                total_tax_ids_abundance[tax_id] = 1
                total_tax_ids_threshold[tax_id] = 1
            results_abundance[tool] = tool_results
            results_read_count[tool] = tool_results
            results_threshold[tool] = tool_results
        else:
            report_tools_abundance[tool] = 1
            report_tools_read_count[tool] = 1  
            tool_results = {}
            tool_results_read_count = {}
            tool_results_threshold = {}
            for read_id in all_transformed_rows[tool]:
                (tax_id, target_rank) = all_transformed_rows[tool][read_id]
                if target_rank != "species":
                    continue
                if read_id in target_reads_sizes:
                    if tax_id not in tool_results:
                        tool_results[tax_id] = 0
                        tool_results_read_count[tax_id] = 0
                    tool_results[tax_id] += target_reads_sizes[read_id]
                    tool_results_read_count[tax_id] += 1
                    if int(target_reads_sizes[read_id]) > threshold:
                        if tax_id not in tool_results_threshold:
                            tool_results_threshold[tax_id] = 0
                        tool_results_threshold[tax_id] += target_reads_sizes[read_id]

            tool_summary = {}
            tool_summary_read_count = {}
            tool_summary_threshold = {}
            tool_sum = 0.0
            tool_sum_read_count = 0.0
            tool_sum_threshold = 0.0

            for tax_id in tool_results:
                if tax_id in genome_sizes:
                    percentage = float(tool_results[tax_id]) / (float(genome_sizes[tax_id]) * 1000000)
                    tool_sum += percentage
            for tax_id in tool_results:
                if tax_id in genome_sizes:
                    abundance = float(tool_results[tax_id]) / (float(genome_sizes[tax_id]) * 1000000)
                    tool_summary[tax_id] = (abundance / tool_sum) * 100.0
                    total_tax_ids_abundance[tax_id] = 1

            for tax_id in tool_results_threshold:
                if tax_id in genome_sizes:
                    percentage = float(tool_results_threshold[tax_id]) / (float(genome_sizes[tax_id]) * 1000000)
                    tool_sum += percentage
            for tax_id in tool_results_threshold:
                if tax_id in genome_sizes:
                    abundance = float(tool_results_threshold[tax_id]) / (float(genome_sizes[tax_id]) * 1000000)
                    tool_summary_threshold[tax_id] = (abundance / tool_sum) * 100.0
                    total_tax_ids_threshold[tax_id] = 1
            
            for tax_id in tool_results_read_count:
                tool_sum_read_count += float(tool_results_read_count[tax_id])

            for tax_id in tool_results_read_count:
                tool_summary_read_count[tax_id] = (tool_results_read_count[tax_id] / tool_sum_read_count) * 100.0
                total_tax_ids_abundance[tax_id] = 1

            results_abundance[tool] = tool_summary
            results_read_count[tool] = tool_summary_read_count
            results_threshold[tool] = tool_summary_threshold

    report_abundance = generate_report(total_tax_ids_abundance, results_abundance)
    abundance_filename = root_reports + "/" + database + "_" + dataset + "_" + target_rank + "_abundance.csv"
    write_report(report_abundance, number_of_nodes_files, taxonomy_names_lists, abundance_filename, report_tools_abundance)

    report_abundance_threshold = generate_report(total_tax_ids_threshold, results_threshold)
    abundance_threashold_filename = root_reports + "/" + database + "_" + dataset + "_" + target_rank + "_abundance_30p.csv"
    write_report(report_abundance_threshold, number_of_nodes_files, taxonomy_names_lists, abundance_threashold_filename, report_tools_abundance)

    report_read_count = generate_report(total_tax_ids_abundance, results_read_count)
    read_count_filename = root_reports + "/" + database + "_" + dataset + "_" + target_rank + "_abundance_read-cnt.csv"
    write_report(report_read_count, number_of_nodes_files, taxonomy_names_lists, read_count_filename, report_tools_read_count)

def main_func(number_of_nodes_files, datasets_path, should_read):
    ranks_lists = {}
    taxonomy_tree_lists = {}
    taxonomy_names_lists = {}

    for nodes_index in range(0, int(number_of_nodes_files)):
        nodes_file = open("nodes" + str(nodes_index) + ".dmp", "r")
        nodes_lines = nodes_file.readlines()
        ranks = {}
        taxonomy_tree = {}
        for line in nodes_lines:
            parts = re.split(r'\t\|\t+', line.strip())
            taxonomy_tree[parts[0].strip()] = parts[1].strip()
            ranks[parts[0].strip()] = parts[2].strip()
        ranks_lists[nodes_index] = ranks
        taxonomy_tree_lists[nodes_index] = taxonomy_tree

    for names_index in range(0, int(number_of_nodes_files)):
        names_file = open("names" + str(names_index) + ".dmp", "r")
        names_lines = names_file.readlines()
        taxonomy_names = {}
        for line in names_lines:
            parts = re.split(r'\|+', line.strip())
            if parts[3].strip() == "scientific name":
                taxonomy_names[parts[0].strip()] = parts[1].strip()
        taxonomy_names_lists[names_index] = taxonomy_names

    root_results = "results"
    root_cleaned_results = "cleaned_results"
    root_reports = "reports"
    genome_sizes_filename = "genome_sizes.txt"

    target_ranks = ["species", "genus"]

    missing_tax_ids_total = {}

    for nodes_index in range(0, int(number_of_nodes_files)):
        missing_tax_ids = {}
        for target_rank in target_ranks:
            missing_file_path = "missing_" + str(target_rank) + str(nodes_index)
            missing_file = open(missing_file_path, "r")
            missing_lines = missing_file.readlines()
            missing = {}
            for line in missing_lines:
                parts = re.split(r'\t+', line.strip())
                missing[parts[0].strip()] = parts[1].strip()
            missing_tax_ids[target_rank] = missing
        missing_tax_ids_total[nodes_index] = missing_tax_ids

    genome_sizes_file = open(genome_sizes_filename, "r")
    genome_sizes_lines = genome_sizes_file.readlines()

    genome_sizes = {}
    for line in genome_sizes_lines:
        parts = re.split(r'\t+', line.strip())
        genome_sizes[parts[0].strip()] = parts[1].strip()

    databases = ["human2", "custom2"]
    datasets = [
        ("01_Mock_100000-bacteria-l1000-q10.fastq", False, 1, False, False),
        ("02_silico-30p-human-70p-bac.fastq", False, 2, False, False),
        ("03_silico-10-bacteria-100k-reads.fasta", False, 3, True, False),
        ("04_silico-3-euka-bac-100k-reads.fastq", False, 4, True, False),
        ("05_human-pathogen.fastq", False, 5, True, False),
        ("06_50-bac-100k.fastq", False, 6, True, False),
        ("08_negative2_10bac_shuffled_human_20k.fasta", False, 8, True, True),
        ("09_zymo_pacbio.fastq", True, 9, False, False),
        ("10_zymo_ont.fastq", True, 10, False, False),
        ("11_SRR11606871_subsampled.fastq", True, 11, False, False),
        ("20_Sample10_ONT_ERR3201941.fastq", True, 20, False, False),
        ("21_Sample20_ONT_ERR3201951.fastq", True, 21, False, False),
        ("22_Sample21_ONT_ERR3201952.fastq", True, 22, False, False),
        ("23_SRR15489011.fastq", True, 23, False, False),
        ("24_SRR15489017.fastq", True, 24, False, False),
        ("25_SRR15489009.fastq", True, 25, False, False)
    ]
    tools = {
        "truth": 0,
        "bracken": 0,
        "centrifuge": 0,
        "clark": 0,
        "clark-s": 0,
        "deSamba": 0,
        "kaiju": 1,
        "kraken": 0,
        "megan": 0,
        "megan-p": 1,
        "metamaps": 0,
        "minimapA": 0,
        "minimapM": 0,
        "ram": 0
    }

    percentage_tools = { "bracken" }
    mapping_tools = { "minimapM", "minimapA", "ram", "deSamba" }

    for (dataset, is_percentage, number_of_dataset, is_kaiju_exception, is_negative) in datasets:
        print("#Dataset: " + str(dataset))
        path_to_dataset = datasets_path + "/" + dataset
        (read_sizes, reads_sizes_kaiju, threshold) = read_read_sizes(dataset, datasets_path)
        for database in databases:
            print("#Database: " + str(database))
            for target_rank in target_ranks:
                print("Target rank: " + str(target_rank))

                all_transformed_rows = {}

                for tool in tools:
                    if should_read:
                        print("#Tool: " + str(tool))
                        is_percentage_tool = False
                        if tool in percentage_tools:
                            is_percentage_tool = False
                        if tool == "truth" and is_percentage:
                            is_percentage_tool = True
                        if tool == "truth" and is_kaiju_exception:
                            transformed_rows = read_cleaned_results(number_of_dataset, database, target_rank, tool+"K", root_cleaned_results, is_percentage_tool)
                            all_transformed_rows[tool+"K"] = transformed_rows
                        transformed_rows = read_cleaned_results(number_of_dataset, database, target_rank, tool, root_cleaned_results, is_percentage_tool)
                        all_transformed_rows[tool] = transformed_rows
                        continue

                    print("#Tool: " + str(tool) + " - start")

                    results_filename = str(root_results) + "/" + str(tool) + "/" + str(database) + "_" + str(number_of_dataset)
                    results_file = open(results_filename, "r") 
                    results_lines = results_file.readlines()

                    parsed_rows = []
                    taxonomy_tree = taxonomy_tree_lists[tools[tool]]
                    ranks = ranks_lists[tools[tool]]

                    parsed_rows_kaiju_exception = []

                    if tool == "centrifuge":
                        parsed_rows = analyse_centrifuge(results_lines)
                    elif tool == "clark":
                        parsed_rows = analyse_clark(results_lines)
                    elif tool == "clark-s":
                        parsed_rows = analyse_clark_s(results_lines)
                    elif tool == "kaiju":
                        parsed_rows = analyse_kaiju(results_lines)
                    elif tool == "kraken":
                        parsed_rows = analyse_kraken(results_lines)
                    elif tool == "megan" or tool == "megan-p" or tool == "megan-d":
                        parsed_rows = analyse_megan(results_lines)
                    elif tool == "metamaps":
                        parsed_rows = analyse_metamaps(results_lines)
                    elif tool == "minimapM" or tool == "ram":
                        parsed_rows = analyse_sam_paf(results_lines, target_rank, taxonomy_tree, ranks, False)
                    elif tool == "minimapA" or tool == "deSamba":
                        parsed_rows = analyse_sam_paf(results_lines, target_rank, taxonomy_tree, ranks, True)
                    elif tool == "truth":
                        parsed_rows = analyse_truth(results_lines, target_rank, taxonomy_tree, ranks, is_percentage)
                        if is_kaiju_exception:
                            parsed_rows_kaiju_exception = analyse_truth_kaiju_exception(results_lines, target_rank, taxonomy_tree, ranks)
                    elif tool == "bracken":
                        parsed_rows = analyse_bracken(results_lines, target_rank, taxonomy_tree, ranks, is_percentage)
                    else:
                        print("Tool " + str(tool) + " not supported")

                    if tool == "truth" and is_kaiju_exception:
                        transformed_rows_kaiju_exception = parsed_rows_kaiju_exception
                        cleaned_filename = root_cleaned_results + "/" + str(tool) + "K_" + str(database) + "_" + str(number_of_dataset) + "_" + str(target_rank) + ".f2"
                        cleaned_outfile = open(cleaned_filename, "w")
                        all_transformed_rows[tool+"K"] = transformed_rows_kaiju_exception
                        for read_id in transformed_rows_kaiju_exception:
                            (tax_id, rank) = transformed_rows_kaiju_exception[read_id]
                            cleaned_outfile.write(read_id.strip() + "\t" + tax_id.strip() + "\t" + rank.strip() + "\n")
                        cleaned_outfile.close()

                    transformed_rows = parsed_rows
                    cleaned_filename = root_cleaned_results + "/" + str(tool) + "_" + str(database) + "_" + str(number_of_dataset) + "_" + str(target_rank) + ".f2"
                    cleaned_outfile = open(cleaned_filename, "w")

                    if tool not in mapping_tools and tool not in percentage_tools and tool != "truth":
                        transformed_rows = transform_results_to_target_rank(parsed_rows, taxonomy_tree, ranks, target_rank)

                    if (tool == "truth" and is_percentage):
                        for tax_id in transformed_rows:
                            percentage = transformed_rows[tax_id]
                            cleaned_outfile.write(tax_id.strip() + "\t" + percentage.strip() + "\n")
                        cleaned_outfile.close()
                    elif tool in percentage_tools:
                        for tax_id in transformed_rows:
                            (percentage, read_count) = transformed_rows[tax_id]
                            cleaned_outfile.write(tax_id.strip() + "\t" + percentage.strip() + "\t" + read_count.strip() + "\n")
                        cleaned_outfile.close()
                    else:
                        for read_id in transformed_rows:
                            (tax_id, rank) = transformed_rows[read_id]
                            cleaned_outfile.write(read_id.strip() + "\t" + tax_id.strip() + "\t" + rank.strip() + "\n")
                        cleaned_outfile.close()

                    all_transformed_rows[tool] = transformed_rows

                generate_read_count_reports(tools, all_transformed_rows, root_reports, database, dataset, target_rank, is_percentage, percentage_tools, is_kaiju_exception, number_of_nodes_files, taxonomy_names_lists, missing_tax_ids_total, is_negative)
                if target_rank == "species":
                    generate_abundance_reports(tools, all_transformed_rows, root_reports, database, dataset, is_percentage, percentage_tools, is_kaiju_exception, number_of_nodes_files, taxonomy_names_lists, read_sizes, reads_sizes_kaiju, threshold, genome_sizes)


if __name__ == '__main__':
    arg_count = 3
    if len(sys.argv) < (arg_count + 1):
        print("Script requires " + str(arg_count) + " arguments")
        exit()
    should_read = True if sys.argv[3] == "Y" else False
    main_func(sys.argv[1], sys.argv[2], should_read)

