#!/usr/bin/env python

import subprocess


IN = open('fastqs.txt', 'r')
for line in IN:
	f = line.strip()
	fprefix = f[:-9]
	command = ('bsub -q long -W 24:00 -n 2 -R "rusage[mem=10000]" -R "span[hosts=1]" -o '
		+ fprefix + '.out -e ' + fprefix + '.err ~/project/Research/Schahram/Schahram-project/split_fastq.py -i '
		+ f +' -l 100000000')
	subprocess.call(command, shell=True)

