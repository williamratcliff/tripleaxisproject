import readncnr3 as readncnr
import numpy as N
import scriptutil as SU
import re
from simple_combine import simple_combine
import copy
import pylab
from findpeak3 import findpeak
from openopt import NLP
import scipy.optimize
import scipy.odr
pi=N.pi


class data_item(object):
    def __init__(self,data,selected=True):
        self.data=data
        self.selected=selected

class Qnode(object):
    def __init__(self,q,th=None,th2th=None,qscans=None,other=None,data=None):
        self.q=q
        self.selected=True
        self.mon0=1.0
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
                self.th.append(data_item(mydata))
                #print 'self.th',self.th
            elif N.abs(mydata.metadata['motor4']['step']-2*mydata.metadata['motor3']['step'])<tol and N.abs(mydata.metadata['motor3']['step'])>tol:
                #print currfile, 'th-2th scan'
                self.th2th.append(data_item(mydata))
            else:
                #print currfile, 'strange scan'
                self.other.append(data_item(mydata))
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
    
    def find_node(self,qcenter):
        
        i=0
        index=False
        for qnode in self.qlist:
                q=qnode.q
                #print 'qcenter',qcenter
                #print 'q',q
                if check_q(q,qcenter):
                    index=i
                    break
                i=i+1
        return index
    
    def condense_node(self,index):
        qnode=self.qlist[index]
        print qnode.q
        #print qnode.th
        
        a3=[]
        counts=[]
        counts_err=[]
        monlist=[]
        for mydataitem in qnode.th:
            mydata=mydataitem.data
            monlist.append(mydata.metadata['count_info']['monitor'])
            counts_err.append(N.array(mydata.data['counts_err']))
            counts.append(N.array(mydata.data['counts']))
            a3.append(N.array(mydata.data['a3']))
        a3_out,counts_out,counts_err_out=simple_combine(a3,counts,counts_err,monlist)
        
        #print a3_out.shape
        #print counts_out.shape
        #print counts_err_out.shape
        qnode.th_condensed={}
        qnode.th_condensed['a3']=a3_out
        qnode.th_condensed['counts']=counts_out
        qnode.th_condensed['counts_err']=counts_err_out
        
        print qnode.th_condensed['counts'].std()
        print qnode.th_condensed['counts'].mean()
        print qnode.th_condensed['counts'].max()
        print qnode.th_condensed['counts'].min()
        if 0:
            pylab.errorbar(a3_out,counts_out,counts_err_out,marker='s',linestyle='None',mfc='black',mec='black',ecolor='black')
            pylab.show()       
        return 
    
    def fit_node(self,index):
        qnode=self.qlist[index]
        print qnode.q
        th=qnode.th_condensed['a3']
        counts=qnode.th_condensed['counts']
        counts_err=qnode.th_condensed['counts_err']
        print qnode.th_condensed['counts'].std()
        print qnode.th_condensed['counts'].mean()
        maxval=qnode.th_condensed['counts'].max()
        minval=qnode.th_condensed['counts'].min()
        diff=qnode.th_condensed['counts'].max()-qnode.th_condensed['counts'].min()\
        -qnode.th_condensed['counts'].mean()
        sig=qnode.th_condensed['counts'].std()
        
        if diff-3*sig>0:
            #the difference between the high and low point and
            #the mean is greater than 3 sigma so we have a signal
            p0=findpeak(th,counts,1)
            print 'p0',p0
            #Area center width Bak
            center=p0[0]
            width=p0[1]
            sigma=width/2/N.sqrt(2*N.log(2))
            Imax=maxval-minval
            area=Imax*(N.sqrt(2*pi)*sigma)
            print 'Imax',Imax
            pin=[area,center,width,0]
            
            
            
            if 0:
                oparam=scipy.odr.Model(gauss)
                mydatao=scipy.odr.RealData(th,counts,sx=None,sy=counts_err)
                myodr = scipy.odr.ODR(mydatao, oparam, beta0=pin)
                myoutput=myodr.run()
                myoutput.pprint()
                pfit=myoutput.beta
            
            if 1:
                
            
            Icalc=gauss(pfit,th)
            
            if 1:
                width_x=N.linspace(p0[0]-p0[1],p0[0]+p0[1],100)
                width_y=N.ones(width_x.shape)*(maxval-minval)/2
                pos_y=N.linspace(minval,maxval,100)
                pos_x=N.ones(pos_y.shape)*p0[0]
                if 1:
                    pylab.errorbar(th,counts,counts_err,marker='s',linestyle='None',mfc='black',mec='black',ecolor='black')
                    pylab.plot(width_x,width_y)
                    pylab.plot(pos_x,pos_y)
                    pylab.plot(th,Icalc)
                    pylab.show()
            
        else:
            #fix center
            #fix width
            print 'no peak'
        
        return


def gauss(p,x):
    #Area center width Bak
    
    x0=p[1]
    width=p[2]
    sigma=width/2/N.sqrt(2*N.log(2))
    area=N.abs(p[0])/N.sqrt(2*pi)/sigma
    background=p[3]
    y=background+area*N.exp(-(0.5*(x-x0)*(x-x0)/sigma/sigma))
    return y

#def chisq(p):
              

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
    Qtree.mon0=mon0
    #flist=flist[0:12]
    for currfile in flist:
        #print 'MAIN READ',currfile
        mydata=mydatareader.readbuffer(currfile)
        mydata.data['counts_err']=N.sqrt(mydata.data['counts'])*mon0/mydata.metadata['count_info']['monitor']
        mydata.data['counts']=N.array(mydata.data['counts'])*mon0/mydata.metadata['count_info']['monitor']
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
    print mydata.data.keys()
    myfilebase='magsc'
    myend='bt9'
    mydirectory=r'c:\ce2rhin8\mar10_2009'
    myfilebaseglob=myfilebase+'*.'+myend
    #print myfilebaseglob
    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    #SU.printr(flist)
    
    qtree=readfiles(flist)
    qtree.condense_node(0)
    qtree.condense_node(1)
    qtree.fit_node(0)
    qtree.fit_node(1)