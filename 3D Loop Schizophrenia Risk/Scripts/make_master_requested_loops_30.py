#!/usr/bin/env python
import argparse
import sys
parser = argparse.ArgumentParser(description='combine master_loops file with requested loop lists')
parser.add_argument('-i', help='input file (merged_loops_30)', type = str, required=True)
args = parser.parse_args()


def get_resolution(line):
	start = int(line.split('\t')[1])
	end = int(line.split('\t')[2])
	r = end - start
	if r != 5000 and r != 10000:
		print 'ERROR: not 5kb or 10kb resolution'
		sys.exit()
	return r



def add_qvals(line, F):
	res = get_resolution(line)
	# Coordinates [chr1, x1, x2, chr2, y1, y2]
	coord = line.split('\t')[:6]
	qvals = []
	for cell in ['Astro2', 'Neu2', 'NPC2']:
		found = False
		REQ = open('../' + cell + '/hiccups/default_MAPQ30_master_loops/requested_list_' + str(res), 'r')
		for reqline in REQ:
			reqcoord = reqline.split('\t')[:6]
			if coord == reqcoord:
				found = True
				qvals.append(float(reqline.split('\t')[17]))
				break
		if not found:
			qvals.append('NA')
		REQ.close()
	if len(qvals) != 3:
		print 'ERROR: something went wrong with grabbing qvals'
		sys.exit()
	else:
		F.write('\t'.join(line.strip().split('\t') + map(str,qvals)) + '\n')

def main():
	index = args.i[-2:]
	OUT = open('master_requested_loops_30_' + index, 'w')
	IN = open(args.i, 'r')

	for i, line in enumerate(IN):
		if i % 100 == 0:
			print 'On row: ' + str(i)
		add_qvals(line, OUT)

	IN.close()
	OUT.close()

if __name__ == '__main__':
	main()