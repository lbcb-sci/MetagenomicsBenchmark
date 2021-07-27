import sys
import re
import os

class Row:
    def __init__(self, is_clasified, read_id, tax_id, length, multiclassification):
        self.is_clasified = is_clasified
        self.read_id = read_id
        self.tax_id = tax_id
        self.length = length
        self.multiclassification = multiclassification

    def to_string(self):
        return str(self.is_clasified) + "\t" + str(self.read_id) + "\t" + str(self.tax_id) + "\t" + str(self.length) + "\t" + str(self.multiclassification) + "\n"

def write_rows(rows, fileout):
    outfile = open(fileout,"w") 
    for row in rows:
        outfile.write(row.to_string())
    outfile.close()

def analyse_metamaps(lines, fileout):
    read_first_line = False
    rows = []
    for line in lines:
        if read_first_line == False:
            read_first_line = True
            continue;
        parts = re.split(r'\t+', line.strip())
        row = Row("C", parts[0].strip(), parts[1].strip(), "", 1.0)
        rows.append(row)
    write_rows(rows, fileout)

def analyse_centrifuge(lines, fileout):
    read_first_line = False
    rows = []
    for line in lines:
        if read_first_line == False:
            read_first_line = True
            continue;
        parts = re.split(r'\t+', line.strip())
        row = Row("C", parts[0].strip(), parts[2].strip(), parts[6].strip(), 1.0)
        rows.append(row)
    write_rows(rows, fileout)

def analyse_clark(lines, fileout):
    read_first_line = False
    rows = []
    for line in lines:
        if read_first_line == False:
            read_first_line = True
            continue;
        parts = re.split(r'\,+', line.strip())
        tax_id = parts[2].strip()
        classified = "C"
        if tax_id == "NA":
            classified = "U"
        row = Row(classified, parts[0].strip(), tax_id, parts[1].strip(), 1.0)
        rows.append(row)
    write_rows(rows, fileout)

def analyse_clark_s(lines, fileout):
    read_first_line = False
    rows = []
    for line in lines:
        if read_first_line == False:
            read_first_line = True
            continue;
        parts = re.split(r'\,+', line.strip())
        tax_id = parts[3].strip()
        classified = "C"
        if tax_id == "NA":
            classified = "U"
        row = Row(classified, parts[0].strip(), tax_id, parts[1].strip(), 1.0)
        rows.append(row)
    write_rows(rows, fileout)


def analyse_kraken(lines, fileout):
    rows = []
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        row = Row(parts[0].strip(), parts[1].strip(), parts[2].strip(), parts[3].strip(), 1.0)
        rows.append(row)
    write_rows(rows, fileout)

def analyse_megan(lines, fileout):
    rows = []
    for line in lines:
        parts = re.split(r'\t+', line.strip())
        row = Row("C", parts[0].strip(), parts[1].strip(), "0", 1.0)
        rows.append(row)
    write_rows(rows, fileout)

def main_func(tool, results_path, fileout):
    results_file = open(results_path, "r") 
    results_lines = results_file.readlines()

    if tool == "kraken":
        analyse_kraken(results_lines, fileout)
    elif tool == "clark":
        analyse_clark(results_lines, fileout)
    elif tool == "clark-s":
        analyse_clark_s(results_lines, fileout)
    elif tool == "centrifuge":
        analyse_centrifuge(results_lines, fileout)
    elif tool == "metamaps":
        analyse_metamaps(results_lines, fileout)
    elif tool == "megan":
        analyse_megan(results_lines, fileout)
    else:
        print("Tool " + str(tool) + " not supported")

if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Script requires 3 arguments")
        exit()
    main_func(sys.argv[1], sys.argv[2], sys.argv[3])

