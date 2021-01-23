import sys
import re
import os

reports_dir = sys.argv[1]

file_write = open("all_tax_in_report","w") 

onlyfiles = [f for f in os.listdir(reports_dir) if os.path.isfile(os.path.join(reports_dir, f))]
for filename in onlyfiles:
	fullpath = reports_dir + "/" + filename
	fullfilename, file_extension = os.path.splitext(fullpath)

	if(file_extension == ".report"):
		file = open(fullpath, "r")
		lines = file.readlines()

		ids = {}

		for line in lines:
			parts = re.split(r'\t+', line.strip())
			ids[parts[0].strip()] = parts[0].strip()

for key in ids:
	file_write.write(str(key.strip()) + "\n")
