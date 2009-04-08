import readncnr3 as readncnr
import numpy as N
import scriptutil as SU
import re
import simple_combine
import copy



class Qnode(object):
    def __init__(self,q,th=None,th2th=None,qscans=None,other=None,data=None):
        self.q=q
        if th==None:
            self.th=[]
        else:
            self.th=th
        if th2th==None:
            self.th2th=[]
        else:
            self.th2th=th2th
        if other==None:
            self.other=[]
        else:
            self.other=other
        if qscans==None:
            self.qscans=[]
        else:
            self.qscans=qscans
        if data!=None:
            self.place_data(data)
        return
    
    def place_data(self,mydata,tol=1e-6):
        if mydata.metadata['file_info']['scantype']=='b':
            #print 'b'
            currfile=mydata.metadata['file_info']['filename']
            if N.abs(mydata.metadata['motor4']['step'])<tol and N.abs(mydata.metadata['motor3']['step'])>tol:
                #print currfile, 'a3 scan'
                self.th.append(mydata)
                #print 'self.th',self.th
            elif N.abs(mydata.metadata['motor4']['step']-2*mydata.metadata['motor3']['step'])<tol and N.abs(mydata.metadata['motor3']['step'])>tol:
                #print currfile, 'th-2th scan'
                self.th2th.append(mydata)
            else:
                #print currfile, 'strange scan'
                self.other.append(mydata)
        return
                
class Qtree(object):
    def __init__(self,qlist=None):
        if qlist==None:
            self.qlist=[]
        else:
            self.qlist=qlist
        return
    def addnode(self,mydata):
        mydata=copy.deepcopy(mydata)
        qcenter=mydata.metadata['q_center']
        qlist=copy.deepcopy(self.qlist)
        if len(qlist)==0:
            #print '0 case'
            newnode=Qnode(qcenter,data=mydata)
            self.qlist.append(newnode)
            #print '0000000000000'
            #print 'qlist len',len(self.qlist)
            #print 'before qlist q',self.qlist[0].q
            #print
        else:
            #print 'else'
            #print 'qlist len',len(self.qlist)
            #print 'before qlist q',self.qlist[0].q
            inlist=False
            for qnode in self.qlist:
                q=qnode.q
                #print 'qcenter',qcenter
                #print 'q',q
                if check_q(q,qcenter):
                    #print qcenter,'in list'
                    qnode.place_data(mydata)
                    #print 'placed'
                    #print qnode.th
                    inlist=True 
            if inlist==False:
                #print 'NOT in list'
                newnode=Qnode(qcenter,data=mydata)
                self.qlist.append(newnode)
        return
            
              

def check_q(q1,q2,tol=1e-6):
    heq=False
    keq=False
    leq=False
    #print 'q1,q2',q1,q2
    if N.abs(q2['h_center']-q1['h_center'])< tol:
        heq=True
    if N.abs(q2['k_center']-q1['k_center'])< tol:
        keq=True
    if N.abs(q2['l_center']-q1['l_center'])< tol:
        leq=True
    #print 'heq',heq,'keq',keq,'leq',leq
    
    return (heq and keq and leq)
   
    

def readfiles(flist,tol=1e-4):
    mydatareader=readncnr.datareader()
    H=[]#N.array([])
    I=[]#N.array([])
    Ierr=[]#N.array([])
    monlist=[]
    count=0
    myfirstdata=mydatareader.readbuffer(flist[0])
    mon0=myfirstdata.metadata['count_info']['monitor']
    print 'mon0',mon0
    qtree=Qtree()
    #flist=flist[0:12]
    for currfile in flist:
        #print 'MAIN READ',currfile
        mydata=mydatareader.readbuffer(currfile)
        qtree.addnode(copy.deepcopy(mydata))
        #print 'readloop'
        #print 'q in loop', qtree.qlist[0].q
        
    for qnode in qtree.qlist:
        print qnode.q['h_center'],qnode.q['k_center'],qnode.q['l_center'],len(qnode.th),qnode.th
    
    #print qtree.qlist
    return qtree

if __name__=='__main__':
    myfilestr=r'C:\Ce2RhIn8\Mar10_2009\magsc035.bt9'
    #myfilestr=r'c:\bifeo3xtal\jan8_2008\9175\fpx53418.bt7'
    #myfilestr=r'c:\13165\13165\data\MagHigh56784.bt7'
    #myfilestr=r'c:\13176\data\CeOFeAs57255.bt7.out'
    mydatareader=readncnr.datareader()
    mydata=mydatareader.readbuffer(myfilestr)
#    print mydata.__dict__
#    print mydata.additional_metadata
#    print mydata.metadata
#    print mydata.metadata['file_info']['scantype']
#    print mydata.metadata['collimations']
#    print mydata.metadata['dspacing']['monochromator_dspacing']
#    print mydata.metadata['dspacing']['analyzer_dspacing']
#    print mydata.metadata['lattice']['a']
#    print mydata.metadata['lattice']['b']
#    print mydata.metadata['lattice']['c']
#    print mydata.metadata['lattice']['alpha']
#    print mydata.metadata['lattice']['beta']
#    print mydata.metadata['lattice']['gamma']
#    print mydata.metadata['motor3']['step']
#    print mydata.metadata['motor4']['step']
#    print mydata.metadata['count_info']['monitor']
#    print mydata.metadata['energy_info']
#    print mydata.metadata['q_center']['h_center']
#    print mydata.metadata['q_center']['k_center']
#    print mydata.metadata['q_center']['l_center']
#    print mydata.metadata['file_info']['filebase']
#    print mydata.metadata['file_info']['filename']
#    print mydata.metadata['file_info']['fileseq_number']
    #print mydata.data.keys()
    myfilebase='magsc'
    myend='bt9'
    mydirectory=r'c:\ce2rhin8\mar10_2009'
    myfilebaseglob=myfilebase+'*.'+myend
    #print myfilebaseglob
    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    #SU.printr(flist)
    
    readfiles(flist)
    
    