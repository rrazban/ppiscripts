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
listi=[6, 7, 8]


def Getaasgen(ver):
	dat = ppisettings.commonseq +ver+ '.dat'
	ophile = open(dat, 'r')
	next(ophile)	#dont have to do awkward iterating
	result = np.zeros((stdline, len(listi)), dtype = "float")
	for l,line in enumerate(ophile):
		for index in listi:
			split = line.split()
			result[l][index - 6] = split[index]
	return result

def CountnFreq(time,resultdict):
	result = np.zeros((len(listi)), dtype = 'float')
	for index in listi:
		for ver in upversion:
			result[index - 6] += float(resultdict[ver][time][index - 6])/float(len(upversion))
	return result

def mp_preentropy(upversion, nprocs):
	def worker(vers, out_q):
		outdict = {}
		for ver in vers:
			outdict[ver] = Getaasgen(ver)
#put countnfreq here		
		out_q.put(outdict)
	
	out_q=Queue()
	chunksize=int(math.ceil(len(upversion)/float(nprocs)))
	procs=[]
	
	for i in range(nprocs):
		p=Process(target=worker,args=(upversion[chunksize*i:chunksize*(i+1)],out_q))
		procs.append(p)
		p.start()

	resultdict={}
	for i in range(nprocs):
		resultdict.update(out_q.get())
	for p in procs:
		p.join()
	return resultdict	

def mp_count(stdline, nprocs, resultdict):
	def worker(vers, out_q):
		outdict={}
		for time in vers:
			outdict[time]=CountnFreq(time,resultdict)
		
		out_q.put(outdict)
	
	out_q=Queue()
	chunksize=int(math.ceil(len(range(stdline))/float(nprocs)))
	procs=[]
	
	for i in range(nprocs):
		p=Process(target=worker,args=(range(stdline)[chunksize*i:chunksize*(i+1)],out_q))
		procs.append(p)
		p.start()

	resultdict={}
	for i in range(nprocs):
		resultdict.update(out_q.get())
	for p in procs:
		p.join()
	return resultdict	

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
	upversion,incomplete,inpresent = ppijobstatus.mp_fail(ppisettings.dirs,nprocs)

	resultdict = mp_preentropy(upversion,nprocs)
	resultdict1=mp_count(stdline,nprocs,resultdict)
	AvgintervalPrinter(resultdict1)

	end=str(datetime.now())
	print 'start time: '+begin
	print 'end time: '+end
