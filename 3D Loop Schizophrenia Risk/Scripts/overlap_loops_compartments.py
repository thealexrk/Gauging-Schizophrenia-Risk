#!/usr/bin/env python

import argparse
import myfunctions as mf
import sys
import gzip



parser = argparse.ArgumentParser(description='Return percentages of A-A, B-B, or A-B loop anchors')
parser.add_argument('-c', help='compartment call file (ex. Neu.553x2_2607.11654.500kb.pc1.bedGraph)',
	type=str, required=True)
parser.add_argument('-t', help='type of cell (ex. Astro, GM, Neu, or NPC)', type=str, required=True)
parser.add_argument('-l', help='cleaned loop file (ex. clean_master_requested_loops)', type=str,
	required=True)
parser.add_argument('-s', help='cell type specific loops', type=bool, default=False)
args=parser.parse_args()

def eigen_to_compartment(e):
	'''
	Return compartment type of 
	eigenvector value
	Args: e float
	Return: c string 
	'''
	if e <= 0 :
		c = 'B'
	elif e > 0 :
		c = 'A'
	else:
		print 'ERROR: undefined eigenvector value'
		sys.exit()
	return c

def get_loop_type(chrom, loop, cf):
	'''
	Get compartment type of loop anchors:
	Args: 
		chrom: string chromosome (ex. chr14)
		loop: tuple of lists with anchor coordinates
		cf: compartment file
	Return: 
		loop_type: string AA, BB, AB, or U (undefined)
	'''
	# Open compartment file
	if cf[-3:] == '.gz':
		C = gzip.open(cf, 'rb')
	else:
		C = open(cf, 'r')

	anchor1_compartment_found = False
	anchor2_compartment_found = False

	for line in C:
		splitline = line.strip().split('\t')
		c_chrom = splitline[0]
		c_start = int(splitline[1])
		c_stop = int(splitline[2])
		eigen = splitline[3]
		if eigen != 'nan':
			eigen = float(eigen)
			eigen = eigen_to_compartment(eigen)
			c_interval = [c_start, c_stop]
			if c_chrom == chrom:
				if mf.contains(c_interval, loop[0]) and not anchor1_compartment_found:
					anchor1_compartment = eigen
					anchor1_compartment_found = True
				if mf.contains(c_interval, loop[1]) and not anchor2_compartment_found:
					anchor2_compartment = eigen
					anchor2_compartment_found = True
	
	if not anchor1_compartment_found or not anchor2_compartment_found:
		loop_type = 'U'
	else:
		loop_type = anchor1_compartment + anchor2_compartment
		if loop_type == 'BA':
			loop_type = 'AB'
	C.close()
	return loop_type
			


def main():

	# Open loop file
	L = open(args.l, 'r')

	# Initialize loop types
	loop_types = {
		'AA': 0,
		'BB': 0,
		'AB': 0,
		'U': 0
	}

	# Parse loop file
	# Skip header
	L.readline()
	for line in L:
		if mf.is_loop_in_cell(line, args.t, args.s, False):
			split_loop_line=line.strip().split('\t')
			loop_chrom = split_loop_line[0]
			loop_start1 = int(split_loop_line[1])
			loop_stop1 = int(split_loop_line[2])
			loop_start2 = int(split_loop_line[4])
			loop_stop2 = int(split_loop_line[5])
			mf.check_start_below_stop(loop_start1, loop_stop1)
			mf.check_start_below_stop(loop_start2, loop_stop2)
			loop1_interval = [loop_start1, loop_stop1]
			loop2_interval = [loop_start2, loop_stop2]
			loop_type = get_loop_type(loop_chrom, (loop1_interval, loop2_interval), args.c)
			loop_types[loop_type] += 1
	if args.s:
		OUT = open('overlap_loops_compartments_' + args.t + '_specific.txt', 'w')
	else:
		OUT = open('overlap_loops_compartments_' + args.t + '.txt', 'w')
	OUT.write('\t'.join(['cell', 'AA', 'BB', 'AB', 'U', 'total']) + '\n')
	OUT.write('\t'.join([args.t, str(loop_types['AA']), str(loop_types['BB']),
		str(loop_types['AB']), str(loop_types['U']), str(sum(loop_types.values()))]) + '\n')
	OUT.close()
	L.close()

		
if __name__ == '__main__':
	main()