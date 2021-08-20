#!/usr/bin/env python
import argparse
import re
import sys

parser=argparse.ArgumentParser(description= 'get distance range for intxns in *-specific_PGC_CLOZUK_intrxns.txt')
parser.add_argument('-f', help='intxns file (*-specific_PGC_CLOZUK_intrxns.txt)', type=str, required=True)
args=parser.parse_args()

def is_cis(a, b):
	asearch  = re.search('chr\d+|chrX|chrY|chrM', a)
	chr1 = asearch.group()
	bsearch = re.search('chr\d+|chrX|chrY|chrM', b)
	chr2 = bsearch.group()
	if chr1 == chr2:
		return True
	else:
		return False

def get_distance(a, b):
	# which bin in first
	astart = int(re.search(':\s*(\d+)', a).group(1))
	bstart = int(re.search(':\s*(\d+)', b).group(1))
	if astart < bstart:
		anchor1 = int(re.search('-\s*(\d+)', a).group(1))
		anchor2 = bstart
	elif bstart < astart:
		anchor1 = int(re.search('-\s*(\d+)', b).group(1))
		anchor2 = astart
	else:
		print 'EQUAL loci!'
		sys.exit()
	d = anchor2 - anchor1
	return d


def main():
	distances = []
	IN=open(args.f, 'r')
	IN.readline()
	for line in IN:
		s = line.split('\t')
		loc1 = s[2]
		loc2 = s[5]
		if is_cis(loc1, loc2):
			d = get_distance(loc1, loc2)
			distances.append(d)
		else:
			print 'Trans interaction!!'
			sys.exit()
	print 'Minimum interaction distance: ' + str(min(distances))
	print 'Maximum interaction distance: ' + str(max(distances))


		



if __name__ == '__main__':
	main()

