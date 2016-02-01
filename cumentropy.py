#!/usr/bin/python

#calculate cumulative sequence entropy and hydrophobicity in parallel

import numpy as np
from multiprocessing import *
from datetime import datetime
import ppijobstatus, ppisettings, translation


aalen = ppisettings.aalen
stdline = ppisettings.stdline 

def readaas(ver):
	aamat=np.empty((stdline,aalen), dtype='float')
	hydromat=np.empty((stdline,aalen), dtype='float')

	dat = ppisettings.commonseq +ver+ '.dat'
	with open(dat) as ophile:
		next(ophile)	#dont have to do awkward iterating
		for l,line in enumerate(ophile):
			split = line.split()
			RNA = split[len(split)-1]
			AAs = translation.rnatoaa(RNA)
			for pos,AA in enumerate(AAs):
				aamat[l][pos] = translation.aatnummap[AA]
				hydromat[l][pos] = translation.aathydromap[AA]
	return aamat,hydromat

def count(currcountmat,aamat,time):
	for pos in range(aalen):
		currcountmat[pos][aamat[time][pos]] += 1
	return currcountmat

def frequency(currcountmat,time):
	return currcountmat/float(time)

def entropy(freq):
	entropy=np.zeros((aalen), dtype='float')
	for pos in range(aalen):
		for num in range(len(translation.aminoacids)):
			if freq[pos][num]==0:
				pass
			else:
				entropy[pos] += -freq[pos][num]*np.log(freq[pos][num])
	return entropy

def mp_cumentropy(versions, nprocs):
	def worker(batch, out_q):
		raw=np.zeros((stdline,2,aalen), dtype='float')
		for ver in batch:
			currcountmat=np.zeros((aalen,len(translation.aminoacids)), dtype='float')
			aamat,hydromat = readaas(ver)
			for	time in range(1,stdline): 
				currcountmat = count(currcountmat,aamat,time)
				freqmat = frequency(currcountmat,time)
				raw[time][0] = np.add(raw[time][0], entropy(freqmat)/float(len(versions)))
				raw[time][1] = np.add(raw[time][1], np.mean(hydromat[:time,:], axis=0)/float(len(versions)))
		out_q.put(raw)
	
	out_q=Queue()
	chunksize=int(np.ceil(len(versions)/float(nprocs)))
	procs=[]
	for i in range(nprocs):
		shape=(aalen)
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
			for l in range(1,stdline):
				for pos in range(aalen):
					wout.write('{0}  {1}  {2}\n'.format(l,pos,data[l][i][pos]))
'''	avgoutputs =['avgcumS.txt', 'avghydro.txt']
	for i,out in enumerate(avgoutputs):
		with open(out, 'w') as wout:
			for l in range(1,stdline):
				wout.write('{0} {1}\n'.format(l,np.mean(data, axis=2)[l][i]))''' #have matlab do


if __name__=='__main__':
	begin=str(datetime.now())
	print "Running cumentropy.py"
	print 'start time: '+begin

	stati = ppijobstatus.mp_jobstatus(ppisettings.dirs,ppisettings.args.nprocs)
	result = mp_cumentropy(stati[0],ppisettings.args.nprocs)
	writer(result, outputs = ['cumS.txt', 'cumhydro.txt'])

	end=str(datetime.now())
	print 'end time: '+end
