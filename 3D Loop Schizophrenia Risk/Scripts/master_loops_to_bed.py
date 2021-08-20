#!/usr/bin/env python
import argparse

parser=argparse.ArgumentParser(description='Convert master_loops file to bed file of all loci') 
parser.add_argument('-i', help= 'master loops file (ex. master_loops)', type=str, required=True)
args=parser.parse_args()

def make_loci_list(FH):
	loci_list = []
	FH.readline()
	for line in FH:
		splitline = line.split()
		loc1 = '\t'.join(splitline[:3])
		loc2 = '\t'.join(splitline[3:6])
		loci_list.append(loc1)
		loci_list.append(loc2)
	return loci_list

def bed_index(mylist):
	bedfmt_list = []
	for loc in mylist:
		splitline = loc.split()
		chrom = splitline[0]
		start = splitline[1]
		end = splitline[2]
		bed_start = str(int(start) - 1)
		bedloc = '\t'.join([chrom, bed_start, end])
		bedfmt_list.append(bedloc)
	return bedfmt_list




def main():

	IN = open(args.i, 'r')
	OUT = open(args.i + '.bed', 'w')
	
	loci_list = make_loci_list(IN)
	unique_loci = list(set(loci_list))
	bed_loci = bed_index(unique_loci)
	for loci in bed_loci:
		OUT.write(loci+'\n')
	IN.close()
	OUT.close()



if __name__ == '__main__':
	main()