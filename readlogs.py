import math

def readlogfile(logfile):
	with open(logfile,'r') as f:
		content = f.readlines()
	livemu = 0.0
	contnum = len(content) -3 -5
	for i in range(5, len(content)-3):
		ticks = content[i].rsplit(",")
		livemu += float(ticks[2])
	ticks = content[-1].rsplit(",")
	if len(ticks) == 2:
		liveticks = int(ticks[0])
		deadticks = int(ticks[1])
		#liveratio = float(liveticks)/float(liveticks+deadticks)
		livetime = liveticks/float(6e8)
		return (livetime,contnum)
	else:
		return (0,0)

def calcrates(dir,numfiles):
	fulltime = 0
	fulltime2 = 0
	total_livetime = 0
	numevts = 0
	for i in range(numfiles):
		infile = str(dir+"/"+str(i)+".log")
		temptimes = readlogfile(infile)
		numevts += temptimes[1]
		total_livetime += temptimes[0]
		fulltime += temptimes[0]/float(numfiles)
		fulltime2 += temptimes[0]*temptimes[0]/float(numfiles)
		#print(str(i)+", Running total"+str(fulltime))
	rate = float(numevts)/total_livetime
	#print("rate "+str(rate))
	ratesig = fulltime2 - fulltime*fulltime
	#print("ratsig: "+str(ratesig))
	ratesig = math.sqrt(ratesig)
	ratesig = rate*ratesig/fulltime
	#print("ratsig: "+str(ratesig))
	outline = "Rate is "+str(rate)+" +/- "+str(ratesig)+" Hz "+" total events is "+str(numevts)+" Total live time is "+str(total_livetime)
	print(outline)
	return rate
