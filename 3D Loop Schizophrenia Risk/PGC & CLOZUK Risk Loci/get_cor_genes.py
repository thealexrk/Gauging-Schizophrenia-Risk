#!/usr/bin/env python
import sys

celltype = sys.argv[1]



def get_COS(g):
	IN = open('COSgenesPositionsMatrix_remove_zero_var.txt', 'r')
	found = False
	for line in IN:
		cg = line.split('\t')[3]
		if cg == g:
			found =True
			break
	IN.close()
	if not found:
		print 'ERROR: no cos row found'
		print g
		return g, False	
	else:
		return line, True

ID_table = open('CLOZUK_PGC_' + celltype + '_ensgID_biomart.txt', 'r')


# OUTFILES
OUT = open('COS_CLOZUK_PGC_' + celltype + '.txt', 'w')

for i, line in enumerate(ID_table):
	if i % 10 == 0:
		print i
	s = line.strip().split('\t')
	gID = s[1]
	result = get_COS(gID)
	if result[1]:
		OUT.write(result[0])

ID_table.close()
OUT.close()