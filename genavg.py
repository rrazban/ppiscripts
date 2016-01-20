#!/usr/bin/python

#extract Pint (0,1) and ppi

import os,sys,glob,re,math
import numpy as np
from multiprocessing import *
import functools as ft
from datetime import datetime
import ppijobstatus, ppisettings 

stdline = ppisettings.stdline
nprocs = ppisettings.args.nprocs

def Getaasgen(ver):
	listi=[6,7,8]
	dir=commonseq+ver
	dat=commondat+ver+'.dat'
	ophile=open(dir+'/'+dat,'r')
	next(ophile)	#dont have to do awkward iterating
	result={}
	for l,line in enumerate(ophile):
		result[l]={}
		for index in listi:
			split=line.split()
			result[l][index]=split[index]
	return result

def CountnFreq(time,resultdict):
	listi=[6,7,8]
	result={}
	for index in listi:
		result[index]=0
		for ver in upversion:
			result[index]+=float(resultdict[ver][time][index])/float(len(upversion))
	return result


def mp_preentropy(upversion, nprocs):
	def worker(vers, out_q):
		outdict={}
		for ver in vers:
			outdict[ver]=Getaasgen(ver)
		
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
		listi=[6,7,8]
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
	output1=open('P0nat.dat','w')
	output2=open('P1nat.dat','w')
	output3=open('ppi.dat','w')
	for r in range(stdline): 
		output1.write('{0}\n'.format(avgentropy[r][6]))
		output2.write('{0}\n'.format(avgentropy[r][7]))
		output3.write('{0}\n'.format(avgentropy[r][8]))




def stupid(count,upver):
	print upver
	print count

if __name__=='__main__':
	begin=str(datetime.now())
	upversion,incomplete,inpresent = ppijobstatus.mp_fail(ppisettings.dirs,nprocs)

	resultdict = mp_preentropy(upversion,nprocs)
	resultdict1=mp_count(stdline,nprocs,resultdict)
	AvgintervalPrinter(resultdict1)

	end=str(datetime.now())
	print 'start time: '+begin
	print 'end time: '+end
