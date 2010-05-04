from ice.event.communication import BroadcastMessageListener
from ice import * 
from ice import Controller
from ice.commands import *
from ice.clientAPI import *

c=ClientAPI.getInstance('me','localhost')

#To use this, simply make an instantiation like:
#g=GetHowLong('qescan'); g.run()


class GetErrs(BroadcastMessageListener):
	def __init__(self):
		c = Controller.getReference()
		comm = c.getCommMgr()
		comm.addMessageListener(self)

	def actionPerformed(self,me):
		data = me.getData()
		print data

g=GetErrs()


class SimpleImmediateCommand(ImmediateCommand):
        def __init__(self,commandstr):
            self.commandstr=commandstr
	def getCommandString(self):
		s = "%s"%(self.commandstr,)
		return s

	def parseSynchronousResponse(self):
		#cid = self.getCommandId()
		#imq = self.getResponseMessageQueue()
		#cmq = mq.getMessagesForAbsCommandId(cid)
		#f = cmq.remove()
		f = self.getResponse()
		#print f 
                #f=str(f.split('\n')[1]).split(':')[1].strip().split()
                #self.soft=f[0]
                #self.hard=f[1]                
                self.result=f



class GetHowlong(ImmediateCommand):
        def __init__(self,scan,overhead=0.0):
            self.scan=scan
	def getCommandString(self):
		s = "Scan Howlong "+self.scan+"-s "+"-o "+self.overhead#overhead is in seconds
		return s

	def parseSynchronousResponse(self):
		#cid = self.getCommandId()
		#imq = self.getResponseMessageQueue()
		#cmq = mq.getMessagesForAbsCommandId(cid)
		#f = cmq.remove()
		f = self.getResponse()
                print f 
                self.f=f

class GetValues(ImmediateCommand):
        def __init__(self,device,overhead=0.0):
            self.device=device
	def getCommandString(self):
		s = "Print "+self.device
		return s

	def parseSynchronousResponse(self):
		#cid = self.getCommandId()
		#imq = self.getResponseMessageQueue()
		#cmq = mq.getMessagesForAbsCommandId(cid)
		#f = cmq.remove()
		f = self.getResponse()
		print f 
                f=str(f.split('\n')[1]).split(':')[1].strip().split()
                self.soft=f[0]
                self.hard=f[1]                
                self.f=f



class Rate(QueuedCommand):
	def getCommandString(self):
		s = "Rate"
		return s

	def parseSynchronousResponse(self):
		cid = self.getCommandId()
		imq = self.getResponseMessageQueue()
		cmq = imq.getMessagesForAbsCommandId(cid)
                #self.cmq=cmq
		f = cmq.remove()
		#f = self.getResponse()
                print "Rate ",f 
                self.f=f
                #Note, the rate command doesn't seem to actually print the bloody rate!!!

class Count(QueuedCommand):
	def __init__(self,device,duration,printflag=True):
		self.device=device
		self.duration=duration
		self.printflag=printflag
	def getCommandString(self):
		s = "Count %s "%(self.device,)
		if self.printflag:
			s=s+'-p '	
		s=s+str(self.duration)
		#print 'string',s
		return s

	def parseSynchronousResponse(self):
		cid = self.getCommandId()
		imq = self.getResponseMessageQueue()
		cmq = imq.getMessagesForAbsCommandId(cid)
		self.cmq=cmq
		self.imq=imq
		#f = cmq.remove()
		##f = self.getResponse()
        #print "Count ",f 
        #self.f=f
                #Note, the rate command doesn't seem to actually print the bloody rate!!!




class GetScans(ImmediateCommand):
	def getCommandString(self):
		s = "Scan List"
		return s

	def parseSynchronousResponse(self):
		#cid = self.getCommandId()
		#imq = self.getResponseMessageQueue()
		#cmq = mq.getMessagesForAbsCommandId(cid)
		#f = cmq.remove()
		f = self.getResponse()
		#print f
                self.f=f
                self.scanlist=[]
                if f:
                   #print 'scanning'
                   scans=f.split('\n')[1:]
                   for scan in scans:
                      #print str(scan).split(':')
                      tmp=str(scan).split(':')
                      #print 'tmp',tmp
                      if tmp[0] is not '' and len(tmp)>1:
                         #print 'len', len(tmp), tmp
                         self.scanlist.append(tmp[1].strip())
                   print self.scanlist
