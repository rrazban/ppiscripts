#!/usr/bin/python

#determine failed jobs, make sure file 1 completes!

import os,sys,glob,re

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--verbose', metavar='True', default='False')
    parser.add_argument('--directory', required=True, metavar='RUNs',  help='dir to look for file')
    parser.add_argument('--delete', metavar='True', default='False', help="delete inpresent directories")
    parser.add_argument('--nprocs', type=int, default=1)
    return parser.parse_args()
args=parse_args()

try:
	os.chdir(args.directory)
except:
	print 'directory not found'
	sys.exit()

dirs=glob.glob('*seqv*')
prenumber=re.findall(r'\d+',dirs[0])
commonseq=dirs[0][:-len(prenumber[len(prenumber)-1])]
commondat='basiclog-'+commonseq

#assume first file completes
dir=commonseq+str(1)
dat=commondat+str(1)+'.dat'
ophile=open(dir+'/'+dat,'r')
tstep=[]
for stdline,R in enumerate(ophile):
	split=R.split()
	tstep.append(split[0])
	aalen=len(split[len(split)-1])/3	#assume divisible by 3

tstep.remove('00000000')
