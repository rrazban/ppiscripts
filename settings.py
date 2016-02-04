#!/usr/bin/python

#set standards from file 1 for ppilib scripts, make sure file 1 sucessfully completes!

import os,sys,glob,re

def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
#    parser.add_argument('--verbose', metavar='True', default=False)  #not available
#    parser.add_argument('--delete', default=False, help="delete inpresent directories")	#problem alleviated with array submission
    parser.add_argument('--nproc', type=int, default=1)
    parser.add_argument('--dir', required=True, metavar='RUNs',  help='dir to look for file')
    return parser.parse_args()
args=parse_args()

try:
	os.chdir(args.dir)
except:
	print 'directory not found'
	sys.exit()

STDdirs = glob.glob('*seqv*')
prenumber = re.findall(r'\d+',STDdirs[0])
STDfile = STDdirs[0][:-(len(prenumber[len(prenumber)-1]) + len('.dat'))]

#assume first file completes
dat = STDfile +str(1)+ '.dat'
tstep=[]
with open(dat) as ophile:
	for STDtimelen,line in enumerate(ophile):
		split = line.split()
		tstep.append(split[0])
		STDseqlen = len(split[len(split)-1])/3	#assume divisible by 3

#tstep.remove('00000000')
