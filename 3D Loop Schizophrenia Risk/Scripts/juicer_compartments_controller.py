#!/usr/bin/env python
import subprocess
import sys

hicfile = sys.argv[1]

chroms = [str(x) for x in range(1,23)] + ['X']

for chrom in chroms:
	command =  'bsub -q short -W 4:00 -R "rusage[mem=20000]" -n 1 -o ' + chrom + '.out ' \
	 	'~/project/Research/Schahram/Schahram-project/juicer_compartments.py ' \
		'-i ' + hicfile + ' -c ' + chrom + ' -r 50000 -g ' \
		'~/project/Research/digest/feature_analysis/gene_density/gene_density_50kb.bedGraph'

	subprocess.call(command, shell=True)

