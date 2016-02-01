#!/usr/bin/python

#calculate cumulative sequence entropy in parallel

import numpy as np
from multiprocessing import *
from datetime import datetime
import ppijobstatus, ppisettings, translation, cumentropy

aalen = ppisettings.aalen
stdline = ppisettings.stdline

def countnfreq(count,aamat,hydromat,versions):
	for time in range(stdline):
		for pos in range(aalen):
			count[0][time][pos][aamat[time][pos]] += 1/float(len(versions))
			count[1][time][pos][0] += hydromat[time][pos]/float(len(versions))
	return count

def mp_preentropy(versions,nprocs):
	def worker(batch, out_q):
		raw=np.zeros((2,stdline,aalen,len(translation.aminoacids)), dtype='float') #notice that only S requires len(translation.aminoacids) dimension
		for ver in batch:
			aamat,hydromat = cumentropy.readaas(ver)
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
	entropy=np.zeros((stdline,2,aalen), dtype='float')	#to be compatible with cumentropy.printer
	for time in range(stdline):
		for pos in range(aalen):
			for num in range(len(translation.aminoacids)):
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

	stati = ppijobstatus.mp_jobstatus(ppisettings.dirs,ppisettings.args.nprocs)
	result = mp_preentropy(stati[0],ppisettings.args.nprocs)
	finalresult = entropy(result)
	cumentropy.writer(finalresult, outputs = ['avgS.txt', 'avghydro.txt'])

	end=str(datetime.now())
	print 'end time: '+end
