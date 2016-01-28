#!/usr/bin/python

#calculate cumulative sequence entropy in parallel

import math
import numpy as np
from multiprocessing import *
from datetime import datetime
import ppijobstatus, ppisettings, translation

aalen = ppisettings.aalen
stdline = ppisettings.stdline
mapaa = translation.mapaa
nprocs = ppisettings.args.nprocs

def initialize():
	counter={}
	for pos in range(aalen):
		counter[pos]={}
		counter[pos]={}
		for key in mapaa:
			counter[pos][mapaa[key]]=0
	return counter

def Getaas(ver):
	shape=(stdline,aalen)
	Bitp=np.empty(shape, dtype='string')

	dat=ppisettings.commonseq +ver+ '.dat'
	ophile=open(dat, 'r')
	next(ophile)	#dont have to do awkward iterating
	for l,line in enumerate(ophile):
		split=line.split()
		RNA=split[len(split)-1]
		AAs=translation.rnatoaa(RNA)
		for pp,AA in enumerate(AAs):
			Bitp[l][pp]=AA
	return Bitp

def mp_preentropy(upversion, nprocs):
	def worker(vers, out_q):
		outdict={}
		for ver in vers:
			outdict[ver]={}
			outdict[ver]=Getaas(ver)
		
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

def CountnFreq(time,resultdict,count):
	for ver in upversion:
		for pos in range(aalen):
			count[pos][resultdict[ver][time][pos]]+=1/float(len(upversion))
	return count


def Entropy(freq):
	entropy={}
	for pos in range(aalen):
		entropy[pos]=0	
		for aa in translation.aminoacids:
			if freq[pos][aa]==0:
				pass
			else:
				entropy[pos]+=-freq[pos][aa]*np.log(freq[pos][aa])
	return entropy


def AvgintervalPrinter(avgentropy):
	output=open('avgS.txt','w')
	output.write('#trange avgentropy\n')
	for r in range(stdline): 
		for pp in range(aalen):
			output.write('{0}  {1}  {2}\n'.format(r,pp,avgentropy[r][pp]))

def PosPrinter(pos,avgentropy):
	output=open('Pos'+str(pos)+'.txt','w')
	for r in range(stdline): 
		output.write('{0}\n'.format(avgentropy[r][pos]))

def mp_count(nprocs, resultdict):
	def worker(vers, out_q):
		outdict={}
		for time in vers:
			count=initialize()
			freq=CountnFreq(time,resultdict,count)
			outdict[time]=Entropy(freq)
		
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


if __name__=='__main__':
	begin=str(datetime.now())
	upversion,incomplete,inpresent=ppijobstatus.mp_fail(ppisettings.dirs,nprocs)

	resultdict=mp_preentropy(upversion,nprocs)
	avgentropy=mp_count(nprocs,resultdict)
	AvgintervalPrinter(avgentropy)
	end=str(datetime.now())
	print 'start time: '+begin
	print 'end time: '+end
