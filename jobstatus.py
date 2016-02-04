#!/usr/bin/python

#determine failed jobs in parallel
# make sure file 1 completes!

import re,math,os,itertools
import numpy as np
from multiprocessing import *
from datetime import datetime
import settings  

stati = ['complete', 'incomplete', 'inpresent']

def mp_jobstatus(files, nprocs):
	nprocs=min(nprocs,2)		 #roughly optimized for 105 jobs, printout 50, mem-per-cpu 100
	def worker(batch, out_q):
		raw=[[] for index in range(len(stati))]
		for file in batch:
			prenumber = re.findall(r'\d+',file)
			number = prenumber[len(prenumber)-1]
			try:
				with open(file) as ofile:
					for nline,line in enumerate(ofile):
						pass
				if nline==settings.STDtimelen:
					raw[0].append(number)
				else:			
					raw[1].append(number)
			except:
				raw[2].append(number)
		out_q.put(raw)

	out_q=Queue()
	chunksize = int(math.ceil(len(files)/float(nprocs)))
	procs=[]
	for i in range(nprocs):
		p=Process(target=worker, args=(files[chunksize*i:chunksize*(i+1)],out_q))
		procs.append(p)
		p.start()

	preresult=[]
	for i in range(nprocs):
		preresult.append(out_q.get())
	for p in procs:
		p.join()
	result=[[] for index in range(len(stati))]
	for i in range(nprocs):
		for index in range(len(stati)):
			result[index].extend(preresult[i][index])
	return result

def printer(result):
	for i,out in enumerate(stati):
		print '{0}: {1} ({2})'.format(out,sorted(result[i]),len(result[i]))

if __name__=='__main__':
	begin=str(datetime.now())
	print 'Running jobstatus.py'
	print 'start time: '+begin

	result = mp_jobstatus(settings.STDdirs,settings.args.nproc)
	printer(result)

	end=str(datetime.now())
	print 'end time: '+end
