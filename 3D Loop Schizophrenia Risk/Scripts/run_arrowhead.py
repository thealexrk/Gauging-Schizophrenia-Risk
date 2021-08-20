#!/usr/bin/env python

import argparse
import subprocess


parser = argparse.ArgumentParser(description='Run Arrowhead TAD calling')
parser.add_argument('-i', help= 'input .hic file', type=str, required=True)
args= parser.parse_args()


def main():
 res = ['10000', '25000', '50000', '100000']
 for r in res:
 	command = ('bsub -q long -W 24:00 -R "rusage[mem=10000]" -o arrowhead/'
 		+ args.i + '_' + r + '.out -e arrowhead/' + args.i + '_' + r + '.err '
 		+ 'java -Xmx2g -jar /home/tb37w/project/Research/Schahram/test/opt/juicer/scripts/juicebox_tools.7.0.jar '
 		+ 'arrowhead --ignore_sparsity -r ' + r + ' ' + args.i + ' arrowhead')
 	subprocess.call(command, shell=True)

if __name__ == '__main__':
	main()

