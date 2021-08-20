#!/usr/bin/env python
import argparse
import subprocess
import numpy as np

parser=argparse.ArgumentParser(description='Quick check for total cis reads in .hic using STRAW')
parser.add_argument('-i', help='input .hic file', type=str, required=True)
args=parser.parse_args()


def main():
	out = []
	for chrom in ['chr' + x for x in map(str, range(1,23)) + ['X', 'Y', 'M']]:
	 	p = subprocess.Popen('/home/tb37w/project/Research/Schahram/test/straw NONE '+ args.i + ' ' + chrom + ' ' + chrom + ' BP 1000000', stdout=subprocess.PIPE, shell=True)
		contacts = [float(x.split('\t')[2]) for x in p.communicate()[0].strip().split('\n')]
		out = out + contacts
	y = np.array(out)
	total_reads = np.sum(y)
	print 'Total cis contacts in ' + args.i + ': ' + str(total_reads)

if __name__ == '__main__':
	main()