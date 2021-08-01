import subprocess
import sys
import re
import os
import shutil
from os import listdir
from os import path
from os.path import isfile, join
from Bio import SeqIO

with open(sys.argv[1], "rU") as handle:
	for record in SeqIO.parse(handle, "fasta"):
		print(record.id)