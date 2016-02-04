#!/usr/bin/python

#extract Pint (0,1) and ppi

import numpy as np
from multiprocessing import *
from datetime import datetime
import jobstatus, settings 

timelen = settings.STDtimelen

#correspond to P0int, P1int, ppi
list = [6, 7, 8]


def mp_genavg(versions, nprocs):
	nprocs=min(nprocs,4)		#roughly optimized for 105 jobs, printout 50, mem-per-cpu 100
	def worker(batch, out_q):
	#get and normalize value
		raw=np.zeros((timelen,len(list)), dtype='float')
		for ver in batch:
			dat = settings.STDfile +ver+ '.dat'
			ophile = open(dat, 'r')
			with open(dat) as ophile:
				next(ophile)	#dont have to do awkward iterating
				for l,line in enumerate(ophile):
					split = line.split()
					for i,index in enumerate(list): 
						raw[l][i] += float(split[index])/len(versions)
		out_q.put(raw)
	
	out_q = Queue()
	chunksize = int(np.ceil(len(versions)/float(nprocs)))
	procs=[]
	for i in range(nprocs):
		p = Process(target = worker, args = (versions[chunksize*i:chunksize*(i+1)], out_q))
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
			result = np.add(result, preresult[i])		
	return result

def writer(data):
	outputs = ['P0nat.txt', 'P1nat.txt', 'ppi.txt']
	for i,out in enumerate(outputs):
		with open(out, 'w') as wout:
			for l in range(timelen): 
				wout.write('{0}\n'.format(data[l][i]))

if __name__=='__main__':
	begin=str(datetime.now())
	print 'Running genavg.py'
	print 'start time: ' + begin

	stati = jobstatus.mp_jobstatus(settings.STDdirs,settings.args.nproc)
	result = mp_genavg(stati[0],settings.args.nproc)
	writer(result)

	end=str(datetime.now())
	print 'end time: '+end
