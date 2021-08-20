#!/usr/bin/env python
import argparse
import subprocess

parser=argparse.ArgumentParser(description = 'Unzip fastq file, split by line number, fix naming, zip up')
parser.add_argument('-i', help='Input gzipped fastq', type=str, required=True)
parser.add_argument('-l', help='Number of lines for split file', type=str, required=True)
args = parser.parse_args()

def main():
	suffix_ext = '.fastq' 
	print "starting"
	p = subprocess.Popen(['gunzip', args.i])
	p.wait()
	print "done unzipping"
	p = subprocess.Popen(['split', '-l', args.l, args.i[:-3], args.i[:-9] + '_'])
	p.wait()
	print "done splitting"
	p = subprocess.Popen('for f in ' + args.i[:-9] + '_*; do mv "$f" "$f.fastq"; done', shell=True)
	p.wait()
	print "done renaming"
	p = subprocess.Popen(['rm', args.i[:-3]])
	p.wait()
	print "done removing orig"
	subprocess.call('gzip ' + args.i[:-9] + '_*', shell=True)
	print "done zipping"
	



if __name__ == '__main__':
	main()