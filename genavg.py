#!/usr/bin/python

#extract Pint (0,1) and ppi

import math
import numpy as np
from multiprocessing import *
from datetime import datetime
import ppijobstatus, ppisettings 

stdline = ppisettings.stdline
nprocs = ppisettings.args.nprocs

#correspond to P0int, P1int, ppi
list= [6, 7, 8]
#currently hard-coded: listi index associated with 0,1,2

def getncount(vers):
	outdict=np.zeros((stdline,len(list)), dtype='float')
	for ver in vers:
		dat = ppisettings.commonseq +ver+ '.dat'
		ophile = open(dat, 'r')
		next(ophile)	#dont have to do awkward iterating
		for l,line in enumerate(ophile):
			split = line.split()
			for index in list: 
				outdict[l][index-6] += float(split[index])/len(upversion)
	return outdict

def AvgintervalPrinter(avgentropy):
	output1=open('P0nat.txt', 'w')
	output2=open('P1nat.txt', 'w')
	output3=open('ppi.txt', 'w')
	for r in range(stdline): 
		output1.write('{0}\n'.format(avgentropy[r][0]))
		output2.write('{0}\n'.format(avgentropy[r][1]))
		output3.write('{0}\n'.format(avgentropy[r][2]))


if __name__=='__main__':
	begin=str(datetime.now())
	print 'Running genavg.py'
	print 'start time: ' + begin

	upversion,incomplete,inpresent = ppijobstatus.mp_fail(ppisettings.dirs,nprocs)

	resultdict = getncount(upversion)
	AvgintervalPrinter(resultdict)

	end=str(datetime.now())
	print 'end time: '+end
