#!/usr/bin/env python

import subprocess


IN = open('fastqs.txt', 'r')
for line in IN:
	f = line.strip()
	fprefix = f[:-9]
	command = ('sbatch --partition=5days --time=48:00:00 --ntasks=2 --mem=10000 --output='
		+ fprefix + '.out --error=' + fprefix + '.err --wrap="~/Schahram-project/split_fastq.py -i '
		+ f +' -l 100000000"')
	subprocess.call(command, shell=True)


