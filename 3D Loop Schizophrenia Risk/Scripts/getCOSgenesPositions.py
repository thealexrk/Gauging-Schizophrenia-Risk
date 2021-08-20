#!/usr/bin/env python
import sys
import subprocess


def get_gene_pos(gfile):
	suffix = gfile[-2:]
	OUT = open('COSgenesPositions_' + suffix, 'w')

	with open(gfile, 'r') as g:
		for gline in g:
			gsplit = gline.split()
			gene = gsplit[0]
			with open('ensGene_hg19.txt', 'r') as f:
				f.readline()
				for line in f:
					splitline = line.split('\t')
					name = splitline[12]
					if name == gene:
						chrom = splitline[2]
						txStart = splitline[4]
						txEnd = splitline[5]
						OUT.write('\t'.join([chrom, txStart, txEnd]) + 
							'\t' + gline)
						break			
	OUT.close()
	
def main():
	get_gene_pos(sys.argv[1])


if __name__ == '__main__':
	main()