#!/usr/bin/python

#calculate cumulative sequence entropy in parallel

import math
import numpy as np
from multiprocessing import *
from datetime import datetime
import ppijobstatus, ppisettings, translation

mapaa = {"UUU":"F", "UUC":"F", "UUA":"L", "UUG":"L",
       "UCU":"S", "UCC":"S", "UCA":"S", "UCG":"S",
       "UAU":"Y", "UAC":"Y", "UAA":"STOP", "UAG":"STOP",
       "UGU":"C", "UGC":"C", "UGA":"W", "UGG":"W",
       "CUU":"L", "CUC":"L", "CUA":"L", "CUG":"L",
       "CCU":"P", "CCC":"P", "CCA":"P", "CCG":"P",
       "CAU":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
       "CGU":"R", "CGC":"R", "CGA":"R", "CGG":"R",
       "AUU":"I", "AUC":"I", "AUA":"I", "AUG":"M",
       "ACU":"T", "ACC":"T", "ACA":"T", "ACG":"T",
       "AAU":"N", "AAC":"N", "AAA":"N", "AAG":"K",
       "AGU":"S", "AGC":"S", "AGA":"S", "AGG":"S",
       "GUU":"V", "GUC":"V", "GUA":"V", "GUG":"V",
       "GCU":"A", "GCC":"A", "GCA":"A", "GCG":"A",
       "GAU":"D", "GAC":"D", "GAA":"E", "GAG":"E",
       "GGU":"G", "GGC":"G", "GGA":"G", "GGG":"G",}
aminoacids=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

aalen = ppisettings.aalen
stdline = ppisettings.stdline
mapaa = translation.mapaa

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

	dir=ppisettings.commonseq+ver
	dat=ppisettings.commondat+ver+'.dat'
	ophile=open(dir+'/'+dat,'r')
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
	output=open('avgS.dat','w')
	output.write('#trange avgentropy\n')
	for r in range(stdline): 
		for pp in range(aalen):
			output.write('{0}  {1}  {2}\n'.format(r,pp,avgentropy[r][pp]))

def PosPrinter(pos,avgentropy):
	output=open('Pos'+str(pos)+'.dat','w')
	for r in range(stdline): 
		output.write('{}\n'.format(avgentropy[r][pos]))

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


def mp_gen(avgentropy,nprocs):
	def worker(vers, out_q):
		outdict={}
		for pos in vers:
			PosPrinter(pos,avgentropy)
		out_q.put(outdict)
	
	out_q=Queue()
	chunksize=int(math.ceil(len(range(aalen))/float(nprocs)))
	procs=[]
	
	for i in range(nprocs):
		p=Process(target=worker,args=(range(aalen)[chunksize*i:chunksize*(i+1)],out_q))
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
	nprocs=8
	upversion,incomplete,inpresent=ppijobstatus.mp_fail(ppisettings.dirs,nprocs)

	resultdict=mp_preentropy(upversion,nprocs)
	avgentropy=mp_count(nprocs,resultdict)
	AvgintervalPrinter(avgentropy)
	mp_gen(avgentropy,nprocs)
	end=str(datetime.now())
	print 'start time: '+begin
	print 'end time: '+end
