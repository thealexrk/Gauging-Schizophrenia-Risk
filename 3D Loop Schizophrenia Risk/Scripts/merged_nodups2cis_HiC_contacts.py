#!/usr/bin/env python
import argparse
import sys

parser=argparse.ArgumentParser(description='Filter merge_nodups.txt removing trans contacts, ' +
	'intra-fragment reads, and reads <= specified MAPQ threshold')
parser.add_argument('-i', help='input merged_nodups.txt', type=str, required=True)
parser.add_argument('-m', help='MAPQ threshold', type=float, default=0.0)
parser.add_argument('-f', help='format (eg. long, medium)', type=str, default='long')
args=parser.parse_args()

def keep_contact(line, threshold):
	''' 
	Return false if trans-contact, intra-fragment read or below MAPQ threshold
	'''
	splitline = line.split()
	frag1 = splitline[3]
	frag2 = splitline[7]
	chr1 = splitline[1]
	chr2 = splitline[5]
	mapq1 = float(splitline[8])
	mapq2 = float(splitline[11])
	if chr1 != chr2:
		return False
	elif frag1 == frag2:
		return False
	elif mapq1 <= threshold or mapq2 <= threshold:
		return False
	else:
		return True

def keep_contact_medium(line, threshold):
	''' 
	Return false if trans-contact, intra-fragment read or below MAPQ threshold
	(medium format)
	'''
	splitline = line.split()
	frag1 = splitline[4]
	frag2 = splitline[8]
	chr1 = splitline[2]
	chr2 = splitline[6]
	mapq1 = float(splitline[9])
	mapq2 = float(splitline[10])
	if chr1 != chr2:
		return False
	elif frag1 == frag2:
		return False
	elif mapq1 <= threshold or mapq2 <= threshold:
		return False
	else:
		return True


def main():
	with open(args.i, 'r') as IN, open('cis_HiC_contacts.txt', 'w') as OUT:
		for i, contact in enumerate(IN):
			if i%100000 == 0:
				print 'On row: ' + str(i)
			if args.f == 'long':
				if keep_contact(contact, args.m):
					OUT.write(contact)
			elif args.f == 'medium':
				# Fix chromosome format
				csplit = contact.split()
				csplit[2] = 'chr' + csplit[2]
				csplit[6] = 'chr' + csplit[6]
				contact = ' '.join(csplit) + '\n'
				if keep_contact_medium(contact, args.m):
					OUT.write(contact)
			else:
				print 'ERROR: unknown format'
				sys.exit()


if __name__ == '__main__':
	main()
