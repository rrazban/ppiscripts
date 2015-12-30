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
aminoacids = translation.aminoacids

def initialize():
	counter={}
	for pos in range(aalen):
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

def Count(count,Bitp,time):
	for pos in range(aalen):
		count[pos][Bitp[time][pos]]+=1
	return count

def Frequency(count,time):
	timez=float(time)+1
	freq={}
	for pos in range(aalen):
		freq[pos]={}
		for aas in aminoacids:
			freq[pos][aas]=count[pos][aas]/float(timez)
	return freq

def Entropy(freq):
	shape=(aalen)
	entropy=np.array([0.0]*shape)
	for pos in range(aalen):
		for aa in aminoacids:
			if freq[pos][aa]==0.0:
				pass
			else:
				entropy[pos]+=-freq[pos][aa]*np.log(freq[pos][aa])
	return entropy

def mp_preentropy(upversion, nprocs):
	def worker(vers, out_q, nproc):
		shape=(aalen)
		outdict={}
		outdict[nproc]={}
		for time in range(stdline):
			outdict[nproc][time]=np.array(shape*[0.0])
			
		for ver in vers:
			count=initialize()
			Bitp=Getaas(ver)
			for	time in range(stdline): 
				count=Count(count,Bitp,time)	
				freq1=Frequency(count,time)
				outdict[nproc][time]=np.add(outdict[nproc][time],Entropy(freq1)/float(len(upversion)))
	
		out_q.put(outdict)
	
	out_q=Queue()
	chunksize=int(math.ceil(len(upversion)/float(nprocs)))
	procs=[]
	
	for i in range(nprocs):
		shape=(aalen)
		p=Process(target=worker,args=(upversion[chunksize*i:chunksize*(i+1)],out_q,i))
		procs.append(p)
		p.start()

	preresultdict={}
	for i in range(nprocs):
		preresultdict.update(out_q.get())
	for p in procs:
		p.join()
	
	resultdict={}
	for time in range(stdline):
		resultdict[time]=np.array(aalen*[0.0])
		for i in range(nprocs):
			resultdict[time]=np.add(resultdict[time],preresultdict[i][time])
		
	return resultdict	


def intervalAverage(entropy):
	avgentropy={}
	for time in range(stdline):
		avgentropy[time]={}		
		for pos in range(aalen):
			avgentropy[time][pos]=0
			for ver in upversion:
				avgentropy[time][pos]+=entropy[ver][time][pos]/float(len(upversion))
	return avgentropy

def AvgintervalPrinter(avgentropy):
    output=open('cumS.dat','w')
    output.write('#trange pos cumentropy\n')
    for r in range(stdline): 
        for pp in range(aalen):
            output.write('{0}  {1}  {2}\n'.format(r,pp,avgentropy[r][pp]))

def PosPrinter(pos,avgentropy):
	output=open('Posc'+str(pos)+'.dat','w')
	for r in range(stdline): 
		output.write('{}\n'.format(avgentropy[r][pos]))



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

	avgentropy=mp_preentropy(upversion,nprocs)
	AvgintervalPrinter(avgentropy)
	mp_gen(avgentropy,nprocs)
	end=str(datetime.now())
	print 'start time: '+begin
	print 'end time: '+end
