import os

def write_hostfile(np):
	f = open("hostfile.txt", 'w')
	n = np/4
	nlist = [n,]*4
	r = np%4
	for i in range(r):
		nlist[i] = nlist[i] + 1
	
	for i in range(4):
		f.write('dgpm-compute-%d.local slots=%d\n' % (i+1, nlist[i]))
	
	f.close()

curdir = os.getcwd()
os.chdir(curdir)

repeat_num = 1
nTotalCores = 48

nlist = [1,2,3,6,12,24,36,48]
nlist = [4,]

n = len(nlist)

res = 10
f = open("%s/result_%dm.txt"%(curdir, res), 'a+')

import time,datetime
if __name__ == "__main__":
	for np in nlist:
		a = [res,]*6
		resultDir = "%s/%dm_n%d" % (curdir, res, np)
		if not os.path.exists(resultDir):
			os.mkdir(resultDir)
		a.append(resultDir)

		write_hostfile(np)
		for i in range(repeat_num):
			t = datetime.datetime.now()
			print t
			s = 'mpiexec -hostfile hostfile.txt ./purposive_sampling /data/NewTestData/tif_%dm/slope_s_%dm.tif,/data/NewTestData/tif_%dm/prof_s_%dm.tif,/data/NewTestData/tif_%dm/wet_s_%dm.tif 2 10 0.01 50 0.7 3 100 %s' % tuple(a)
			print s
			result = os.popen(s)
			print result.read()
			t1=datetime.datetime.now()
			print t1
			dt=t1-t
			#print result.read()
			print np, dt
			f.write("%d %s" % (np, dt))
		f.write('\n')
		f.flush()
	f.close()
