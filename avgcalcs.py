#!/usr/bin/python

#calculate average sequence entropy and hydrophobicity in parallel
import numpy as np
from multiprocessing import *
from datetime import datetime
import settings,jobstatus,maps,cumcalcs

seqlen = settings.STDseqlen
timelen = settings.STDtimelen

def countnfreq(count,aamat,hydromat,versions):
	for time in range(timelen):
		for pos in range(seqlen):
			count[0][time][pos][aamat[time][pos]] += 1/float(len(versions))
			count[1][time][pos][0] += hydromat[time][pos]/float(len(versions))
	return count

def mp_preentropy(versions,nprocs):
	nprocs=min(nprocs,32)		 #roughly optimized for 105 jobs, printout 50, mem-per-cpu 100
	def worker(batch, out_q):
		raw=np.zeros((2,timelen,seqlen,len(maps.aminoacids)), dtype='float') #notice that only S requires len(maps.aminoacids) dimension
		for ver in batch:
			aamat,hydromat = cumcalcs.readseq(ver)
			raw = countnfreq(raw,aamat,hydromat,versions)
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


def entropy(freq):
	entropy=np.zeros((timelen,2,seqlen), dtype='float')	#shape chosen to be compatible with cumcalcs.printer
	for time in range(timelen):
		for pos in range(seqlen):
			for num in range(len(maps.aminoacids)):
				if freq[0][time][pos][num]==0:
					pass
				else:
					entropy[time][0][pos] += -freq[0][time][pos][num]*np.log(freq[0][time][pos][num])
				entropy[time][1][pos] = freq[1][time][pos][0]
	return entropy

if __name__=='__main__':
	begin=str(datetime.now())
	print 'Running avgentropy'
	print 'start time: '+begin

	stati = jobstatus.mp_jobstatus(settings.STDdirs,settings.args.nproc)
	result = mp_preentropy(stati[0],settings.args.nproc)
	finalresult = entropy(result)
	cumcalcs.writer(finalresult, outputs = ['avgS.txt', 'avghydro.txt'])

	end=str(datetime.now())
	print 'end time: '+end
