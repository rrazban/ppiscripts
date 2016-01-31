#!/usr/bin/python

#extract Pint (0,1) and ppi

import math
import numpy as np
from multiprocessing import *
from datetime import datetime
import ppijobstatus, ppisettings 

stdline = ppisettings.stdline
nprocs = ppisettings.args.nprocs

#correspond to P0int, P1int, ppi
list = [6, 7, 8]
#currently hard-coded for diff of 1 between ea element 
#list index associated with 0,1,2


def mp_preentropy(upversion, nprocs):
	def worker(vers, out_q):
	#get and normalize value
		out=np.zeros((stdline,len(list)), dtype='float')
		for ver in vers:
			dat = ppisettings.commonseq +ver+ '.dat'
			ophile = open(dat, 'r')
			with open(dat) as ophile:
				next(ophile)	#dont have to do awkward iterating
				for l,line in enumerate(ophile):
					split = line.split()
					for index in list: 
						out[l][index-min(list)] += float(split[index])/len(upversion)
		out_q.put(out)
	
	out_q = Queue()
	chunksize = int(math.ceil(len(upversion)/float(nprocs)))
	procs=[]
	for i in range(nprocs):
		p = Process(target = worker, args = (upversion[chunksize*i:chunksize*(i+1)], out_q))
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

def printer(data):
	outputs = ['P0nat.txt', 'P1nat.txt', 'ppi.txt']
	for index,out in enumerate(outputs):
		with open(out, 'w') as wout:
			for l in range(stdline): 
				wout.write('{0}\n'.format(data[l][index]))

if __name__=='__main__':
	begin=str(datetime.now())
	print 'Running genavg.py'
	print 'start time: ' + begin

	upversion,incomplete,inpresent = ppijobstatus.mp_fail(ppisettings.dirs,nprocs)

	result = mp_preentropy(upversion,nprocs)
	printer(result)

	end=str(datetime.now())
	print 'end time: '+end
