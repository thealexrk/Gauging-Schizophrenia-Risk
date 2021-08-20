#!/usr/bin/env python
import myfunctions as mf

import argparse

parser=argparse.ArgumentParser(description='overlap loops and eQTL')
parser.add_argument('-i', help='loop file (ex. clean_master_requested_loops)', type=str, required=True)
parser.add_argument('-o', help='out file', type=str, required=True)
parser.add_argument('-e', help= 'eQTL file (ex. cis_eQTL_PGC_SCZ2_filter_FDR_1e-20.tsv)', type=str, required=True)
args=parser.parse_args()



def main():
	loop_file = args.i
	LFH = open(loop_file, 'r')
	eQTL_file = args.e
	outfile = args.o
	OUT = open(outfile, 'w')
	# Output file header
	#OUT.write(LFH.readline().strip() + '\tGENE\tGENE_CHR\tGENE_BP1\tGENE_BP2\tSNP\tSNP_CHR\tSNP_POS\n')
	
	for i, line in enumerate(LFH):
		EFH = open(eQTL_file, 'r')
		EFH.readline()
		gene_names = []
		gene_chroms = []
		gene_starts = []
		gene_stops = []
		snp_names = []
		snp_chroms = []
		snp_poss = []
		if i%20 == 0:
			print 'On row: ' + str(i)
		splitline = line.split('\t')
		loop_chrom = splitline[0]
		loop_start1 = int(splitline[1])
		loop_stop1 = int(splitline[2])
		loop_start2 = int(splitline[4])
		loop_stop2 = int(splitline[5])
		mf.check_start_below_stop(loop_start1, loop_stop1)
		mf.check_start_below_stop(loop_start2, loop_stop2)
		loop1_interval = [loop_start1, loop_stop1]
		loop2_interval = [loop_start2, loop_stop2]
		for eline in EFH:
			esplitline = eline.split('\t')
			gene_chrom = 'chr' + esplitline[4]
			gene_start = int(esplitline[5])
			gene_stop = int(esplitline[6])
			mf.check_start_below_stop(gene_start, gene_stop)
			gene_interval = [gene_start, gene_stop]
			snp_chrom = 'chr' + esplitline[8]
			snp_pos = int(esplitline[9])

			# Check if loop and eQTL overlap
			if gene_chrom == snp_chrom == loop_chrom:
				# Holy shit here comes an ugly comparison... 
				if ((mf.overlap(loop1_interval, gene_interval) 
					and (loop_start2 <= snp_pos <= loop_stop2)) 
					or (mf.overlap(loop2_interval, gene_interval)
					and (loop_start1 <= snp_pos <= loop_stop1))):
						gene_names.append(esplitline[2])
						gene_chroms.append(gene_chrom)
						gene_starts.append(str(gene_start))
						gene_stops.append(str(gene_stop))
						snp_names.append(esplitline[7])
						snp_chroms.append(snp_chrom)
						snp_poss.append(str(snp_pos))
		
		OUT.write(line.strip() + '\t' + '|'.join(gene_names) +
			'\t' + '|'.join(gene_chroms) + '\t' + '|'.join(gene_starts) +
			'\t' + '|'.join(gene_stops) + '\t' + '|'.join(snp_names) +
			'\t' + '|'.join(snp_chroms) + '\t' + '|'.join(snp_poss) + '\n')
		EFH.close()

if __name__ == '__main__':
	main()
