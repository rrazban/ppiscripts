#!/usr/bin/python

#calculate cumulative sequence entropy and hydrophobicity in parallel

import numpy as np
import gzip, random
from multiprocessing import *
from datetime import datetime
import settings, jobstatus, maps

seqlen = settings.STDseqlen
timelen = settings.STDtimelen

def readseq(ver):
	aamat=np.empty((timelen,seqlen), dtype='float')
	hydromat=np.empty((timelen,seqlen), dtype='float')

	dat = settings.STDfile +str(ver)+ '.dat'
	with open(dat) as ophile:
		for l,line in enumerate(ophile):
			split = line.split()
			RNA = split[len(split)-1]
			AAs = maps.rnatoaa(RNA)
			for pos,AA in enumerate(AAs):
				aamat[l][pos] = maps.aatnum[AA]
				hydromat[l][pos] = maps.aathydro[AA]
	return aamat,hydromat

def count(currcountmat,aamat,time):
	for pos in range(seqlen):
		currcountmat[pos][aamat[time][pos]] += 1
	return currcountmat

def frequency(currcountmat,time):
	return currcountmat/float(time+1)

def entropy(freq):
	entropy=np.zeros((seqlen), dtype='float')
	for pos in range(seqlen):
		for num in range(len(maps.aminoacids)):
			if freq[pos][num]==0:
				pass
			else:
				entropy[pos] += -freq[pos][num]*np.log(freq[pos][num])
	return entropy

def mp_cumentropy(versions, nprocs):
	nprocs=min(nprocs,32)		 #roughly optimized for 105 jobs, printout 50, mem-per-cpu 100
	def worker(batch, out_q):
		raw=np.zeros((timelen,2,seqlen), dtype='float')
		for ver in batch:
			currcountmat=np.zeros((seqlen,len(maps.aminoacids)), dtype='float')
			aamat,hydromat = readseq(ver)
			for	time in range(timelen): 
				currcountmat = count(currcountmat,aamat,time)
				freqmat = frequency(currcountmat,time)
				raw[time][0] = np.add(raw[time][0], entropy(freqmat)/float(len(versions)))
				raw[time][1] = np.add(raw[time][1], np.mean(hydromat[:(time+1),:], axis=0)/float(len(versions)))
		out_q.put(raw)
	
	out_q=Queue()
	chunksize=int(np.ceil(len(versions)/float(nprocs)))
	procs=[]
	for i in range(nprocs):
		p=Process(target=worker,args=(versions[chunksize*i:chunksize*(i+1)],out_q))
		procs.append(p)
		p.start()

	preresult=[]
	for i in range(nprocs):
		preresult.append(out_q.get())
	for p in procs:
		p.join()
	for i in range(nprocs):
		if i==0:
			result = preresult[i]
		else:						
			result = np.add(result,preresult[i])
	return result

def writer(data, outputs):
	for i,out in enumerate(outputs):
		with open(out, 'w') as wout:
			for l in range(timelen):
				for pos in range(seqlen):
					wout.write('{0}  {1}  {2}\n'.format(l,pos,data[l][i][pos]))
'''	avgoutputs =['avgcumS.txt', 'avghydro.txt']
	for i,out in enumerate(avgoutputs):
		with open(out, 'w') as wout:
			for l in range(1,timelen):
				wout.write('{0} {1}\n'.format(l,np.mean(data, axis=2)[l][i]))''' #have matlab do


if __name__=='__main__':
	begin=str(datetime.now())
	print "Running cumcalcs.py"
	print 'start time: '+begin

	stati = jobstatus.mp_jobstatus(settings.STDdirs,settings.args.nproc)
#	stati[0] = random.sample(stati[0], 100)
#	print stati[0]
	result = mp_cumentropy(stati[0],settings.args.nproc)
	writer(result, outputs = ['cumS.txt', 'cumhydro.txt'])

	end=str(datetime.now())
	print 'end time: '+end
