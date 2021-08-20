#!/usr/bin/env python
import argparse
import sys

parser=argparse.ArgumentParser(description='Filter HiCPro validPairs file removing trans ' +
	'contacts and format for Juicer pre')
parser.add_argument('-i', help='input *_allValidPairs', type=str, required=True)
args=parser.parse_args()


def is_chrom(c1, c2):
	'''
	Check strings are chromosomes
	'''
	if c1[:3] != 'chr' or c2[:3] != 'chr':
		return False
	else:
		return True

def keep_contact(line):
	''' 
	Return false if trans-contact
	'''
	splitline = line.split('\t')
	chr1 = splitline[1]
	chr2 = splitline[4]
	if not is_chrom(chr1, chr2):
		print 'ERROR selecting chromosome'
		sys.exit()
	if chr1 != chr2:
		return False
	else:
		return True
def pre_format(line):
	'''
	Format contact line in short format for input to Juicer pre command
	'''
	splitline = line.split('\t')
	strd_map = {'+': '0', '-': '1'}
	strd1 = splitline[3]
	strd2 = splitline[6]
	new_strd1 = strd_map[strd1]
	new_strd2 = strd_map[strd2]
	chr1 = splitline[1]
	chr2 = splitline[4]
	if not is_chrom(chr1, chr2):
		print 'ERROR selecting chromosome'
		sys.exit()
	pos1 = str(int(splitline[2]))
	pos2 = str(int(splitline[5]))
	frag1 = '0'
	frag2 = '1'
	p = [new_strd1, chr1, pos1, frag1, new_strd2, chr2, pos2, frag2]
	pf = '\t'.join(p) + '\n'
	return pf




def main():
	with open(args.i, 'r') as IN, open('cis_HiC_contacts.txt', 'w') as OUT:
		for i, contact in enumerate(IN):
			if i%100000 == 0:
				print 'On row: ' + str(i)
			if keep_contact(contact):
				pre_contact = pre_format(contact)
				OUT.write(pre_contact)



if __name__ == '__main__':
	main()
