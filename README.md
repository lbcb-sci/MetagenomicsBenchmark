
# Metagenomics benchmark

### Legend

The commands executed to obtain results from the tools are listed in the `all_commands.txt`

`nucleotide_database_filelist.txt` contains list of directories containing genome files in the nucleotide database.

`genome_sizes.txt` - average lengths of genomes obtained from the NCBI web page.
`missing_<rank><database_type>` - files containing <rank> (species/genus) that are in the datasets but are not present in the database for nucleotide (database_type=0) and protein (database_type=1) databases.

The results are calculated for each dataset where every dataset has its order number.
The numbers for each dataset:

- 1 - 01_Mock_100000-bacteria-l1000-q10.fastq
- 2 - 02_silico-30p-human-70p-bac.fastq
- 3 - 03_silico-10-bacteria-100k-reads.fasta
- 4 - 04_silico-3-euka-bac-100k-reads.fastq
- 5 - 05_human-pathogen.fastq
- 6 - 06_50-bac-100k.fastq
- 8 - 08_negative2_10bac_shuffled_human_20k.fasta
- 9 - 09_zymo_pacbio.fastq
- 10 - 10_zymo_ont.fastq
- 11 - 11_SRR11606871_subsampled.fastq
- 20 - 20_Sample10_ONT_ERR3201941.fastq
- 21 - 21_Sample20_ONT_ERR3201951.fastq
- 22 - 22_Sample21_ONT_ERR3201952.fastq
- 23 - 23_SRR15489011.fastq
- 24 - 24_SRR15489017.fastq
- 25 - 25_SRR15489009.fastq

Database have these names:
 - custom - database built from refseq Bacteria and Archaea genomes
 - human - same database as custom + human genome

File types:
 - .f2 files - cleaned classification results obtained after processing tools results.
 - each row in the .f2 file is in format: `<read_id>	<tax_id>	<percentage>	<level>`
 	- <level> - taxon level: species, genus, family...
 - report files - report file, every row contains number of reads classified to a taxon: `<tax_id>	<species_name>	<true_number>	<kraken_number>	<centrifuge_number>	<clark_number>	<metamaps_number>	<megan_number>	...`

### Supporting files:

Download the supporting files from the following link: 
https://zenodo.org/deposit/7198915

The supporting files contain:

1. results - results of tools without minimap2 and ram results because they were too big to upload. We can provide them on request. Truth results included.
2. cleaned_results - results of the classification after the parsing of tool results, stored in .f2 format Cleaned results are calculated for human and custom database and for species and genus level. The format of file names is `<tool>_<database>_<dataset_number>_<rank>.f2`.

3. nodes.dmp and names.dmp taxonomy files
4. reports - Containig report results. Files are named in format: `<database>_<dataset_number>_<dataset_name>_<rank>_<report_type>.csv`
	- <format type> - can be: 
		- `read_count`, 
		- `read_count_TN` (true negatives), 
		- `read_count_TP` (true positives), 
		- `abundance`, 
		- `abundance_30p` (abundance with only 30% longest reads), 
		- `abundance_read-cnt` 
5. sequences - dataset sequences

When downloading supporting files and placing them in the root directory the benchmark scrypt can be run with:

`python benchmark.py 2 <path_to_sequences> <should_read_cleaned_results>` 

`<should_read_cleaned_results>` params determines if the reports will be generated from cleaned results all from original results. Value `Y` will generate reports from the cleaned results, every other value will generate it from original results. In order to generate reports from original results either ram and minimap results need to be provided or those tools can be commented out in the scrypt.   
