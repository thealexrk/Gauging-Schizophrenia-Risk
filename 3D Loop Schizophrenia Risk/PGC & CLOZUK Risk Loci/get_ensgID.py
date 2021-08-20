#!/usr/bin/env python
import sys

def get_ID(g):
	found = False
	e = 0
	IN = open('/home/Tyler/Research/Schahram/genes/ensembl/mart_export_hg38.txt', 'r')
	for line in IN:
		ensGene = line.strip().split('\t')[1]
		if ensGene == g:
			e = line.split('\t')[0]
			found = True 
			break
	IN.close()

	if not found:
		print 'ERROR: no id found'
		print g
		
	return e



CLZ = open('CLOZUK_lists_for_RNA.txt', 'r')
CLZ.readline()

# OUTFILES
CLZ_NPC = open('CLOZUK_NPC_ensgID_biomart.txt', 'w')
CLZ_Neu = open('CLOZUK_Neu_ensgID_biomart.txt', 'w')
CLZ_Glia = open('CLOZUK_Glia_ensgID_biomart.txt', 'w')

# Missing genes from ENSEMBL
missing_NPC = open('missing_NPC_ensgID_biomart.txt', 'w')
missing_Neu = open('missing_Neu_ensgID_biomart.txt', 'w')
missing_Glia = open('missing_Glia_ensgID_biomart.txt', 'w')

for i, line in enumerate(CLZ):
	# NPC
	if i % 10 == 0:
		print i
	s = line.strip('\n').split('\t')
	NPC_g = s[0]
	if NPC_g != '':
		gID = get_ID(NPC_g)
		if gID == 0:
			missing_NPC.write(NPC_g + '\n')
		else:
			CLZ_NPC.write(NPC_g + '\t' + gID + '\n')
	Neu_g = s[1]
	if Neu_g != '':
		gID = get_ID(Neu_g)
		if gID == 0:
			missing_Neu.write(Neu_g + '\n')
		else:
			CLZ_Neu.write(Neu_g + '\t' + gID + '\n')
	Glia_g = s[2]
	if Glia_g != '':
		gID = get_ID(Glia_g)
		if gID == 0:
			missing_Glia.write(Glia_g + '\n')
		else:
			CLZ_Glia.write(Glia_g + '\t' + gID + '\n')

CLZ.close()

# OUTFILES
CLZ_NPC.close()
CLZ_Neu.close()
CLZ_Glia.close()

# Missing genes from ENSEMBL
missing_NPC.close()
missing_Neu.close()
missing_Glia.close()


