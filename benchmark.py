import parse_tool_output
import analize_tool_output
import analize_results
import os
import sys

if __name__ == '__main__':
	# tools = ["kraken", "centrifuge", "clark", "metamaps", "megan"]
	tools = ["kraken"]
	databases = ["default", "custom", "human"]
	technologies = ["ont", "pb"]
	root_parsed = "parsed_results/"
	root_cleaned = "tax_cleaned_results/"
	root_reports = "reports/"

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
				for num in range(1, 5):
					for tech in technologies:
						if num > 2 and tech == "ont":
							continue

						results_filename = "results/" + str(tool) + "/" + str(database) + "_" + tech + "_" + str(num)

						if mode != "only_cleaning":
							parsed_dir = root_parsed + str(tool)
							if os.path.isdir(parsed_dir) == False:
								os.mkdir(parsed_dir)
							parsed_filename = parsed_dir + "/" + str(database) + "_" + tech + "_" + str(num) + ".stat"
							print("Parsing: " + str(tool) + " - " + str(database) + " - " + tech + " - " + str(num))
							parse_tool_output.main_func(tool, results_filename, parsed_filename)

						if mode != "only_parsing":
							cleaned_dir = root_cleaned + str(tool)
							if os.path.isdir(cleaned_dir) == False:
								os.mkdir(cleaned_dir)
							cleaned_filename = cleaned_dir + "/" + str(database) + "_" + tech + "_" + str(num) + ".f2"
							print("Cleaning: " + str(tool) + " - " + str(database)+ " - " + tech + " _ " + str(num))
							analize_tool_output.main_func(tool, parsed_filename, nodes_lines, cleaned_filename)

	if mode != "clean_and_parse" and mode != "only_cleaning" and mode != "only_parsing":
		for database in databases:
			for num in range(1, 5):
				for tech in technologies:
					if num > 2 and tech == "ont":
						continue
					dataset = tech + "_" + str(num)
					print("Analysing: " + str(database) + " - " + dataset)
					analize_results.main_func(root_cleaned, root_reports, dataset, database, names_lines)
