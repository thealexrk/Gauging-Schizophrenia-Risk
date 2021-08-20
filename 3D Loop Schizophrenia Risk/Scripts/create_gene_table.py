#!/usr/bin/env python
import argparse
import sys

parser=argparse.ArgumentParser(description='add gene position to FPKM table')
parser.add_argument('-i', help='FPKM table', type=str, required=True)
args = parser.parse_args()

def main():
	OUT = open('gene_table_' + args.i[-2:], 'w')
	with open(args.i, 'r') as f:
		for line in f:
			gene = line.split('\t')[1]
			if gene[:4] != 'ENSG':
				print 'ERROR: Non-valid gene id'
				sys.exit()
			with open('ensGene_hg19.txt', 'r') as e:
				for eline in e:
					esplitline = eline.split('\t')
					egene = esplitline[12]
					if gene == egene:
						chrom = esplitline[2]
						start = esplitline[4]
						end = esplitline[5]
						OUT.write('\t'.join([chrom, start, end])+'\t'+line)
						break

if __name__ == '__main__':
	main()

