#!/usr/bin/env python
import argparse
import os
import sys
import subprocess

parser=argparse.ArgumentParser(description='HICCUPS submitter for GHPCC')
parser.add_argument('-d', help='work directory under /juicer/work/ (ex. NPC)', type=str, required=True)
parser.add_argument('-m', help='MAPQ > value (ex. 0 or 30)', type=str, required=True)
parser.add_argument('-l', help='loop list with parent directory (ex. master_loops_Schahram/master_loops_30)', type=str, default=False)
args=parser.parse_args()

def make_sh(wkdir, data_dir, mapq, loop_list, hic):
	if loop_list:
		sh_file = data_dir+'_MAPQ'+mapq+'_master_loops_run_hiccups.sh'
		SH = open(sh_file, 'w')
		SH.write(
		'#!/bin/bash\n' +
		'#BSUB -q gpu\n' + 
		'#BSUB -W 2400\n' +
		'#BSUB -R "rusage[mem=20000,ngpus_excl_p=1] span[hosts=1]"\n' +
		'#BSUB -o '+wkdir+data_dir+'/hiccups/default_MAPQ'+mapq+'_master_loops/run_hiccups.out\n' +
		'#BSUB -e '+wkdir+data_dir+'/hiccups/default_MAPQ'+mapq+'_master_loops/run_hiccups.err\n' + 
		'module unload java\n' +
		'module load jdk/1.8.0_31\n' +
		'module load gcc/4.7.4\n' +
		'module load cuda/7.0.28\n' +
		'echo "modules loaded successfully"\n' + 
		'bash juicebox.sh hiccups -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000 '+wkdir+data_dir+'/aligned/'+hic+
		' '+wkdir+data_dir+'/hiccups/default_MAPQ'+mapq+'_master_loops '+wkdir+args.l+'\n'
		)
	else:
		sh_file = data_dir+'_MAPQ'+mapq+'_run_hiccups.sh'
		SH = open(sh_file, 'w')
		SH.write(
		'#!/bin/bash\n' +
		'#BSUB -q gpu\n' + 
		'#BSUB -W 2400\n' +
		'#BSUB -R "rusage[mem=20000,ngpus_excl_p=1] span[hosts=1]"\n' +
		'#BSUB -o '+wkdir+data_dir+'/hiccups/default_MAPQ'+mapq+'/run_hiccups.out\n' +
		'#BSUB -e '+wkdir+data_dir+'/hiccups/default_MAPQ'+mapq+'/run_hiccups.err\n' + 
		'module unload java\n' +
		'module load jdk/1.8.0_31\n' +
		'module load gcc/4.7.4\n' +
		'module load cuda/7.0.28\n' +
		'echo "modules loaded successfully"\n' + 
		'bash juicebox.sh hiccups -r 5000,10000 -f 0.1,0.1 -p 4,2 -i 7,5 -d 20000,20000 '+wkdir+data_dir+'/aligned/'+hic+
		' '+wkdir+data_dir+'/hiccups/default_MAPQ'+mapq+'\n'
		)
	SH.close()
	return sh_file

def check_hicfile(mapq):
	if mapq == '0':
		hicfile = 'inter.hic'
	elif mapq == '30':
		hicfile = 'inter_30.hic'
	else:
		'ERROR: wrong MAPQ value'
		sys.exit()
	return hicfile

def main():
	wkdir = '/home/tb37w/project/Research/Schahram/test/opt/juicer/work/'
	if args.l:
		if not os.path.exists(wkdir+args.d+'/hiccups/default_MAPQ'+args.m+'_master_loops'):
			os.makedirs(wkdir+args.d+'/hiccups/default_MAPQ'+args.m+'_master_loops')
	else:	
		if not os.path.exists(wkdir+args.d+'/hiccups/default_MAPQ'+args.m):
			os.makedirs(wkdir+args.d+'/hiccups/default_MAPQ'+args.m)
	hicfile = check_hicfile(args.m)
	shell_file = make_sh(wkdir, args.d, args.m, args.l, hicfile)
	subprocess.call('bsub < ' + shell_file, shell=True)

	
if __name__ == '__main__':
	main()