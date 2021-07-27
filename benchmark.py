import parse_tool_output
import analize_tool_output
import analize_results
import analize_true_positives
import os
import sys

if __name__ == '__main__':
	tools = ["kraken", "centrifuge", "clark", "metamaps", "megan", "clark-s"]
	databases = ["human2"]
	root_parsed = "parsed_results4/"
	root_cleaned = "tax_cleaned_results4/"
	root_reports = "reports4/"

	start = 3
	number = 4

	mode = sys.argv[1]

	if os.path.isdir(root_parsed) == False:
		os.mkdir(root_parsed)

	if os.path.isdir(root_cleaned) == False:
		os.mkdir(root_cleaned)

	if os.path.isdir(root_reports) == False:
		os.mkdir(root_reports)

	nodes_file = open("nodes.dmp", "r")
	nodes_lines = nodes_file.readlines()

	names_file = open("names.dmp", "r")
	names_lines = names_file.readlines()
	
	if mode != "only_reports":
		for tool in tools:
			for database in databases:
				for num in range(start, number):
					results_filename = "results4/" + str(tool) + "/" + str(database) + "_" + str(num)

					parsed_dir = root_parsed + str(tool)
					if os.path.isdir(parsed_dir) == False:
						os.mkdir(parsed_dir)
					parsed_filename = parsed_dir + "/" + str(database) + "_" + str(num) + ".stat"
					if mode != "only_cleaning":
						print("Parsing: " + str(tool) + " - " + str(database) + " - " + str(num))
						parse_tool_output.main_func(tool, results_filename, parsed_filename)

					if mode != "only_parsing":
						cleaned_dir = root_cleaned + str(tool)
						if os.path.isdir(cleaned_dir) == False:
							os.mkdir(cleaned_dir)
						cleaned_filename = cleaned_dir + "/" + str(database) + "_" + str(num) + ".f2"
						print("Cleaning: " + str(tool) + " - " + str(database)+ " - " + str(num))
						analize_tool_output.main_func(tool, parsed_filename, nodes_lines, cleaned_filename)

	if mode != "clean_and_parse" and mode != "only_cleaning" and mode != "only_parsing" and mode != "vibrio":
		for database in databases:
			for num in range(start, number):
				dataset = str(num)
				print("Analysing: " + str(database) + " - " + dataset)
				analize_results.main_func(root_cleaned, root_reports, dataset, database, names_lines)
				if num < 9:
					analize_true_positives.main_func(root_cleaned, root_reports, dataset, database, names_lines)
