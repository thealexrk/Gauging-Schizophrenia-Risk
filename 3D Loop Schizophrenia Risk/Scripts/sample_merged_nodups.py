#!/usr/bin/env python
import argparse
import random

parser=argparse.ArgumentParser(description='Randomly sample from cis_HiC_contacts.txt')
parser.add_argument('-i', help='input cis_HiC_contacts.txt', type=str, required=True)
parser.add_argument('-k', help='sample size', type=int, required=True)
args = parser.parse_args()

def total_lines(file):
	'''
	Count total lines of file
	'''
	with open(file, 'r') as f:
		for i, line in enumerate(f):
			pass
	return i + 1

def random_sample(file, s, k):
	'''
	Randomly sample k lines from file with s total lines
	without replacement and write sample to file
	'''
	OUT = open(file[:-4] + '_sample_' + str(k) + '.txt', 'w')
	sample_list = sorted(random.sample(xrange(s), k), reverse=True)
	sample = sample_list.pop()
	with open(file, 'r') as f:
		for i, line in enumerate(f):
			if i == sample:
				OUT.write(line)
				if len(sample_list) > 0:
					sample = sample_list.pop()
				else:
					break

def main():
	# Count total number of contacts
	total_contacts = total_lines(args.i)

	# Randomly sample and write to file
	random_sample(args.i, total_contacts, args.k)


if __name__ == '__main__':
	main()