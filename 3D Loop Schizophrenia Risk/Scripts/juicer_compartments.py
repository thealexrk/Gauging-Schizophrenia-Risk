#!/usr/bin/env python

import argparse
import subprocess
import numpy as np
import sys

parser=argparse.ArgumentParser(description='Wrapper for juicer_tools.jar eigenvector to compute eigenvector ' +
	'for chromosome and assign (+/-) sign based on gene density')
parser.add_argument('-i', help='input .hic file', type=str, required=True)
parser.add_argument('-c', help='chromosome (ex. 1)', type=str, required=True)
parser.add_argument('-r', help='resolution to call compartments (ex. 50000)', type=str, default='50000')
parser.add_argument('-g', help='gene density file which must have same resolution (ex. gene_density_50kb.bedGraph)', 
	type=str, required=True)
args=parser.parse_args()

def get_chrom_gene_density(f, c):
	'''
	Get gene density for chromosome from file
	
	Args:
		f: string gene density file
		c: string chromosome
	Returns:
		g: list of gene density for chromosome
		s: list of start bins
		e: list of end bins
	'''
	g = []
	s = []
	e = []
	IN = open(f, 'r')
	for line in IN:
		splitline = line.strip().split('\t')
		chrom = splitline[0]
		if chrom == c:
			g.append(splitline[3])
			s.append(splitline[1])
			e.append(splitline[2])
	IN.close()
	return g, s, e

def check_eigen_sign(g, e):
	'''
	Equivalent of Bryan Lajoie's detectActiveCompartment()
	function from: 
	https://github.com/dekkerlab/cworld-dekker/blob/master/scripts/python/matrix2EigenVectors.py
	
	Args:
		g: list of gene density for chromosome
		e: list of eigenvector values for chromosome
	Returns:
		ne: new list of eigenvector values with sign flipped if necessary
	'''

	posSum = 0
	negSum = 0

	for i in range(len(g)):
		eigen = e[i]
		if eigen == 'NaN':
			continue
		eigen = float(eigen)
		nGenes = float(g[i])
		if eigen > 0:
			posSum += (nGenes*abs(eigen))
		if eigen < 0:
			negSum += (nGenes*abs(eigen))
	if negSum > posSum:
		np_e = np.array(map(float, e))
		ne = map(str, np_e * -1)
	else:
		ne = e
	return ne


def main():
	# Calculate eigenvector for chrom
	command = 'java -jar ~/project/Research/Schahram/test/juicer_tools.1.8.9_jcuda.0.8.jar ' \
		'eigenvector -p KR ' + args.i + ' ' + args.c + ' BP ' + args.r
	p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	eigen_stdout = p.communicate()[0]
	eigen = eigen_stdout.split()

	# Get gene density for chrom
	gene_density, starts, ends = get_chrom_gene_density(args.g, 'chr' + args.c)
	if len(gene_density) != len(eigen):
		print 'ERROR: unequal row numbers'
		sys.exit()

	# Make sure eigenvector + and - signs match gene density (A = + and B = -)
	new_eigen = check_eigen_sign(gene_density, eigen)

	# Write eigenvector to file
	OUT = open(args.i[:-4] + '_chr' + args.c +'_eigen1.bedGraph', 'w')
	for i in range(len(new_eigen)):
		if new_eigen[i] == 'NaN':
			new_eigen[i] = 'nan'
		OUT.write('chr' + args.c +'\t' + starts[i] +'\t'+ ends[i] + '\t' + new_eigen[i] + '\n')
	OUT.close()

if __name__ == '__main__':
	main()
