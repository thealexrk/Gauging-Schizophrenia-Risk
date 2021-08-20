#!/usr/bin/env python
import sys


def overlap(a, b):
        '''
        Return boolean for whether interval b
        overlaps interval a
        '''
        if (a[0] <= b[0] <= a[1]) or (b[0] <= a[0] <= b[1]):
                return True
        else:
                return False
                
def contains(a, b):
	'''
	Return boolean for whether interval a contains
	interval b
	'''
	if (a[0] <= b[0] <= a[1]) and (a[0] <= b[1] <= a[1]):
		return True
	else:
		return False


def check_start_below_stop(start, stop):
	if start >= stop:
		print 'ERROR: start >= stop'
		sys.exit()


def is_loop_in_cell(line, cell, specific, filtr):
	col_id = {'Astro': 21,
				'GM': 22,
				'Neu': 23,
				'NPC': 24}
	splitline = line.rstrip().split('\t')
	qval = float(splitline[col_id[cell]])
	if specific:
		other_cells = col_id.keys()
		other_cells.remove(cell)
		qval1 = float(splitline[col_id[other_cells[0]]])
		qval2 = float(splitline[col_id[other_cells[1]]])
		qval3 = float(splitline[col_id[other_cells[2]]])
		if filtr:
			reads = float(splitline[7])
			if qval < -1 and qval1 >= -1 and qval2 >= -1 and qval3 >= -1 and reads > 50:
				print line.strip()
				return True
			else:
				return False
		else:

			if qval < -1 and qval1 >= -1 and qval2 >= -1 and qval3 >= -1:
				return True
			else:
				return False
	else:
		if qval < -1:
			return True
		else:
			return False

