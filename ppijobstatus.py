#!/usr/bin/python

#determine failed jobs in parallel
# make sure file 1 completes!

import re,math,os,subprocess
from multiprocessing import *
from datetime import datetime
import ppisettings  

def mp_fail(multiple, nprocs):
	def worker(multiple, out_q):
		outdict={}
		outdict['complete']=[]
		outdict['incomplete']=[]
		outdict['inpresent']=[]
		for single in multiple:
			prenumber=re.findall(r'\d+',single)
			number=prenumber[len(prenumber)-1]
			try:
				ofile=open(single+'/'+ppisettings.commondat+number+'.dat','r')		
				for nline,R in enumerate(ofile):
					pass
				if nline==ppisettings.stdline:
					outdict['complete'].append(number)
				else:			
					outdict['incomplete'].append(number)
			except:
				outdict['inpresent'].append(number)
		out_q.put(outdict)


	out_q=Queue()
	chunksize=int(math.ceil(len(multiple)/float(nprocs)))
	procs=[]
	
	for i in range(nprocs):
		p=Process(target=worker,args=(multiple[chunksize*i:chunksize*(i+1)],out_q))
		procs.append(p)
		p.start()

	resultdict=[]
	for i in range(nprocs):
		resultdict.append(out_q.get())
	for p in procs:
		p.join()

	complete=[y for nproc in resultdict for y in nproc['complete']]
	incomplete=[y for nproc in resultdict for y in nproc['incomplete']]
	inpresent=[y for nproc in resultdict for y in nproc['inpresent']]
	return complete,incomplete,inpresent


if __name__=='__main__':
	begin=str(datetime.now())
	nprocs=8
	complete,incomplete,inpresent=mp_fail(ppisettings.dirs,nprocs)
	print 'complete: {0} ({1})'.format(sorted(complete),len(complete))
	print 'incomplete: {0} ({1})'.format(sorted(incomplete),len(incomplete))
	print 'inpresent: {0} ({1})'.format(sorted(inpresent),len(inpresent))
	if ppisettings.args.delete:
		for num in inpresent:
			os.system("rm -r " + ppisettings.commonseq + num)	
		#	subprocess.call(["rm -r ", ppisettings.commonseq, num])
	end=str(datetime.now())
	print 'start time: '+begin
	print 'end time: '+end
