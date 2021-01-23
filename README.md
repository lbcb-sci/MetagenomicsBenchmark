
**Metagenomics benchmark.**


All the results and analysis files are in the following format database_dataset.extension

database:
 - default - database downloaded with or provided by the tools
 - custom - database built from refseq Bacteria and Archaea genomes
 - human - same database as custom + human genome

dataset:
 - ont_x <x = [1:2]> - Mock and 30p ont datasets
 - pb_x <x = [1:4]> - 4 pacbio datasets

extension:
 - .stat - file extension for classification results obtained after processing of classificators results
 	- each row in the .stat file is in format: <is_classified>	<read_id>	<tax_id>	<length>	<percentage>
 	- <is_classified> - C for classified read_id, U for unclassified
 	- <read_id> - id of the classified read
 	- <tax_id> - tax id of the clade to which the read is classified
 	- <length> - length of the read if accessible
 	- <percentage> - some tools classify reads to more than one tax is, if so those reads are classified in <percentage> to certain tax ids. 
 - .f2 - file extension for cleaned classification results obtained after processing .stat files. All rows contain only classified reads.
 	- each row in the .f2 file is in format: <read_id>	<tax_id>	<percentage>	<level>
 	- <level> - taxon level: species, genus, family...
 - .report - report file, every row contains number of reads classified to a taxon: <tax_id>	<species_name>	<true_number>	<kraken_number>	<centrifuge_number>	<clark_number>	<metamaps_number>	<megan_number>



**Contents:**

1. truth - ground truth files in .f2 format
2. results - results of the classifications for evey tool: kraken, clark, metamaps, megan and centrifuge
3. parsed_results - results of the classification after the parsing of tools outputs, stored in .stat format
4. tax_cleaned_results - results of the classitifation after the parsing of .stat files, stored in .f2 format
5. reports - resulting reports stored in .report format

6. nodes.dmp and names.dmp taxonomy files

7. summary_b.txt, summary_a.txt, database_summary.txt
 - summary files of the databases, copied from the refseq Archaea (a) and Bacteria (b) databases. Two a and b summaries are merged in database_summary.txt.

8. genome_sizes_a.txt, genome_sizes_b.txt, database_genome_sizes.txt
 - average lengths of genomes calculated from the database sequences for Arcahea (a), Bacteria (b) and both (database_genome_sizes).

9. all_tax_in_report
 - The collection of tax ids in all reports

10. time-memory
 - files with execution time and memory consumption

**Scripts:**
 
 - parse_tool_output.py - parses classification output files from tools and produces read_id to tax_id mappings in .stat file
 	- arguments: <tool> <results_path> <fileout>
 		- <tool> name of the tool, kraken, clark, centrifuge, metamaps or megan
 		- <results_path> path to the classification results file
 		- <fileout> path to the resulting .stat file

- analize_tool_output.py - takes .stat file and produces .f2 file
	- arguments: <tool> <results_path> <nodes_file> <fileout>
		- <tool> name of the tool
		- <results_path> .stat file
		- <nodes_file> path to the nodes.dmp file
		- <fileout> path to the resulting .f2 file

- analize_results.py - analizes .f2 results from all the tools for one dataset and one database and generates .report file
	- arguments: <root_cleaned> <root_reports> <dataset> <database> <names_file>
		- <root_cleaned> - path to the root folder that contains .f2 files
		- <root_reports> - path to the root folder where the reports are stored
		- <dataset> - name of the dataset
		- <database> - name of the database
		- <names_file> - path to the names.dmp file

- benchmark.py - script that executes whole analysis pipeline, for every tool it executes parse_tool_output, then analize_tool_output and then for every dataset and database it generates reports with analize_results
	- arguments: <mode> there are 4 modes in which this script can be run:
		- only_parsing - only parse_tool_output is executed for every tool, database and dataset, .stat files are generated
		- only_cleaning - only analize_tool_output is executed for evey tool, database and dataset, .f2 files are generated, but the precondition is that .stat files have been generated in advance
		- clean_and_parse - parse_tool_output and analize_tool_output are executed without generating reports
		- only_reports - only .report files are generated, but the precondition is that .f2 files were generated in advance
		- all - all the steps are executed and resulting .report files are generated for evey tool, dataset and database

- collect_all_tax_ids.py - helper script that collects all unique tax_ids from the root directory of stored reports
	- arguments: <reports_dir> - root directory of the reports 
	- tax ids are stored in <all_tax_in_report> file which is initially in this repository

- calculate_average_genomes_lengths.py - helper script that calculates average length of every genome in the database, from the seqeunces in the database
	- genome_sizes_a.txt and genome_sizes_b.txt files were generated with this script using refseq Arcaea and Bacteria databases
	- database_genome_sizes.txt containes both genome_sizes_a.txt and genome_sizes_b.txt
	- arguments: <database_path> <names_path> <filename_out>
		- <database_path> - path to the database containing genome sequences all in the root level directory
		- <names_path> - path to the names.dmp file
		- <filename_out> - path to the output file where the average genome lengths are stored

- fetch_missing_tax_id_lengths.py - helper script that fetches length of any tax id that is not in the database (if that for some reason happens with some tool). The script looks into all_tax_in_report and database_genome_sizes and for every tax id in all_tax_in_report that is not found in database_genome_sizes it fetches the genome sequences for that tax id and calculates average genome length for that tax id
	- arguments: <genome_sizes_file> <tax_ids_file> <summary_file> <outfile>
		- <genome_sizes_file> - file containing genome lengths, for this example database_genome_sizes.txt
		- <tax_ids_file> - file containing all the tax ids in the report files, for this example all_tax_in_report
		- <summary_file> - database summary file, for this example database_summary.txt
		- <outfile> - path to the resulting output file where the lengths of the missing tax ids are stored

- analize_abundances.py - script that prints the .report file containing analysis of abundances







