from __future__ import with_statement
import numpy as N
import matplotlib
from matplotlib import pylab
import paramiko
from contextlib import contextmanager
import socket,sys,os



def fgetl(myfile,returnline=None):
    lineStr=myfile.readline()
    if returnline !=None:
        try:
            returnline[0]=lineStr.rstrip()
        except:
            pass
    strippedLine=lineStr.rstrip()
    #tokenized=strippedLine.split()
    
    return strippedLine

class Data(object):
    def __init__(self):
        pass
    
    
def reader2(myfilestr):
    myfile=open(myfilestr)
    data=Data()
    sR='RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR'
    sA='AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA'
    sI='IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'
    sF='FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF'
    sS='SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS'

    while 1:
        linestr=fgetl(myfile)
        if not(linestr):
            break
        #print linestr
        #get name of file
        if linestr.find(sR)>=0:
            linestr=fgetl(myfile)
            toks=linestr.split()
            data.file_number=int(toks[0])
        if linestr.find(sA)>=0:
            linestr=fgetl(myfile)
            linestr=fgetl(myfile)
            linestr=fgetl(myfile) #we are now on the comment line
            linestr=fgetl(myfile)
            toks=linestr.split()
            data.scan_type=toks[-1]
            
        if linestr.find('        H (Hmin)        K (Kmin)        L (Lmin)             phi             chi')>=0:
            fgetl(myfile);fgetl(myfile)
            fgetl(myfile);fgetl(myfile)
            fgetl(myfile);fgetl(myfile)
            fgetl(myfile);fgetl(myfile);fgetl(myfile)
            hkl=fgetl(myfile).split()
            #Q=sscanf(hkl,'%g');
            data.Hmin=float(hkl[0])
            data.Kmin=float(hkl[1])
            data.Lmin=float(hkl[2])
            data.phi=float(hkl[3])
            data.chi=float(hkl[4])
            hkl=fgetl(myfile).split()
            data.omega=float(hkl[0])
            data.tth=float(hkl[1])
            data.psi=float(hkl[2])
            data.ubmatrix=N.zeros((3,3),'Float64')
            data.ubmatrix[0,0]=float(hkl[3])
            data.ubmatrix[0,1]=float(hkl[4])            
            hkl=fgetl(myfile).split()
            data.ubmatrix[0,2]=float(hkl[0])
            data.ubmatrix[1,0]=float(hkl[1])
            data.ubmatrix[1,1]=float(hkl[2])
            data.ubmatrix[1,2]=float(hkl[3])
            data.ubmatrix[2,0]=float(hkl[4])
            hkl=fgetl(myfile).split()
            data.ubmatrix[2,1]=float(hkl[0])            
            data.ubmatrix[2,2]=float(hkl[1])
            data.wavelength=float(hkl[2])
            data.dmonochromator=float(hkl[3])
            data.danalyzer=float(hkl[4])
    #H (Hmin)        K (Kmin)        L (Lmin)             phi             chi
    #omega  2theta (gamma)             psi         ub(1,1)         ub(1,2)
    #ub(1,3)         ub(2,1)         ub(2,2)         ub(2,3)         ub(3,1)
    #ub(3,2)         ub(3,3)      wavelength  dmonochromator       danalyser
    #energy            Hmax            Kmax            Lmax          DeltaH
    #DeltaK          DeltaL     Deltaenergy         Ki (Kf)       Ddetector
    #2the_ana         ome_ana        samptabl       omega_mon         (spare)
    #scan start       scan step      scan width          preset    add.bkg.step
    #add.bkg.width  add.bkg.preset  couplingfactor         (spare)         (spare)
    #    Temp-s.pt      Temp-Regul     Temp-sample       Voltmeter       Mag.field
            hkl=fgetl(myfile).split()
            data.energy=float(hkl[0])
            data.Hmax=float(hkl[1])
            data.Kmax=float(hkl[2])
            data.Lmax=float(hkl[3])
            data.DeltaH=float(hkl[4])
            hkl=fgetl(myfile).split()
            data.DeltaK=float(hkl[0])
            data.DeltaL=float(hkl[1])
            data.DeltaEnergy=float(hkl[2])
            data.KiKf=float(hkl[3])
            data.Ddetector=float(hkl[4])
            hkl=fgetl(myfile).split()
            data.analyzer_tth=float(hkl[0])
            data.analyzer_omega=float(hkl[1])
            data.sample_table=float(hkl[2])
            data.monochromator_omega=float(hkl[3])
            hkl=fgetl(myfile).split()
            data.scan_start=float(hkl[0])
            data.scan_step=float(hkl[1])
            data.scan_width=float(hkl[2])
    #scan start       scan step      scan width          preset    add.bkg.step
    #add.bkg.width  add.bkg.preset  couplingfactor         (spare)         (spare)
    #    Temp-s.pt      Temp-Regul     Temp-sample       Voltmeter       Mag.field
            data.preset=float(hkl[3])
            data.additional_background_step=float(hkl[4])
            hkl=fgetl(myfile).split()
            data.coupling_factor=float(hkl[2]) # th 2 th?
            hkl=fgetl(myfile).split()
            data.temp_setpoint=float(hkl[0])
            data.temp_regulator=float(hkl[1])
            data.temp_sample=float(hkl[1])
            #let's start reading in real data!
            data.frames=[]
            data.angle1=[]
            data.angle2=[]
            data.monitors=[]
            data.times=[]
            #i=0
            while 1:
                linestr=fgetl(myfile)
                if not(linestr):
                    break
                if linestr.find(sS)>=0:
                    hkl=fgetl(myfile).split()
                    currpt=int(hkl[0])-1
                    data.frames.append([])
                    data.planned_pts=int(hkl[2])
                    
                if linestr.find(sF)>=0:
                        hkl=fgetl(myfile).split()#throw away
                        hkl=fgetl(myfile).split() #throw away
                        hkl=fgetl(myfile).split() 
                        data.times.append(float(hkl[0]))
                        data.monitors.append(float(hkl[1]))
                        data.angle1.append(float(hkl[3])/1000)
                        try:
                            data.angle2.append(float(hkl[4])/1000)
                        except:
                            pass
                if linestr.find(sI)>=0:      
                        i=0
                        hkl=fgetl(myfile).split()  #we know that the detector=32x32 has 1024 points
                        npts=int(hkl[0])
                        flag=True
                        while flag:
                            if i==npts:
                                break
                            hkl=fgetl(myfile).split()
                            for val in hkl:
                                data.frames[-1].append(float(val))
                                i=i+1
                        data.frames[-1]=N.array(data.frames[-1],'float64').reshape((32,32))
                            
        
    return data
                            
                            
                        
            
            
            
 
            
            
        
            
            
        

def reader(myfilestr):
    myfile=open(myfilestr,'r')
    data=Data()
    while 1:
        linestr=fgetl(myfile)
        print linestr
        if not(linestr):
            break
        sst='IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII';
        sss='FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF';
        if linestr.find(sst)>=0: 
            fgetl(myfile);fgetl(myfile);fgetl(myfile);fgetl(myfile);fgetl(myfile);
            hkl=fgetl(myfile);
            Q=hkl.split()
            data.numpoints=float(Q[5])
        if linestr.find('        H (Hmin)        K (Kmin)        L (Lmin)             phi             chi')>=0:
            fgetl(myfile);fgetl(myfile)
            fgetl(myfile);fgetl(myfile)
            fgetl(myfile);fgetl(myfile)
            fgetl(myfile);fgetl(myfile);fgetl(myfile)
            hkl=fgetl(myfile).split()
            #Q=sscanf(hkl,'%g');
            data.Hmin=float(hkl[0])
            data.Kmin=float(hkl[1])
            data.Lmin=float(hkl[2])
            fgetl(myfile);fgetl(myfile);fgetl(myfile)
            hkl=fgetl(myfile).split()
            data.Hmax=float(hkl[1])   
            data.Kmax=float(hkl[2])
            data.Lmax=float(hkl[3])
            fgetl(myfile);fgetl(myfile);fgetl(myfile);fgetl(myfile)
            Temp=fgetl(myfile).split()
            #TTemp=sscanf(Temp,'%g');
            data.Tset=float(Temp[0])
            data.Treg=float(Temp[1])
            data.Tsamp=float(Temp[2])
            #fprintf('Q = ( %3g , %3g , %3g )\n',Qh,Qk,Ql);
            #fprintf('Sample Temperature = %3.2g Kelvin \n',data.Tsamp);
            arr=[]
            fgetl(myfile); fgetl(myfile); fgetl(myfile); fgetl(myfile)
            while 1:
                intensityl=fgetl(myfile)         
                if intensityl!='':
                    print intensityl
                    intensity=intensityl.split()
                    for intens in intensity:
                        arr.append(float(intens))
                else:
                    truelen=len(arr)/8*8
                    data.nm=truelen/8
                    data.intensity=N.array(arr[0:truelen]).reshape((truelen/8,8))
                    Hstep=(data.Hmax-data.Hmin)/(data.numpoints-1);
                    Kstep=(data.Kmax-data.Kmin)/(data.numpoints-1);
                    Lstep=(data.Lmax-data.Lmin)/(data.numpoints-1);
                    Hmax=data.Hmin+(data.nm-1)*Hstep
                    Kmax=data.Kmin+(data.nm-1)*Kstep
                    Lmax=data.Lmin+(data.nm-1)*Lstep
                    data.H=N.linspace(data.Hmin,Hmax,truelen/8)
                    data.K=N.linspace(data.Kmin,Kmax,truelen/8)
                    data.L=N.linspace(data.Lmin,Lmax,truelen/8)
                    return data
            break
    
    #temp=fscanf(fid,'%c',[inf]);data_end=length(temp);
    #intensity=sscanf(temp(1,1:data_end-1),'%g');[n m]=size(intensity);nm=n*m/8
    #myfile.close()
    #data.nm=nm;
    
    
    
@contextmanager
def create_ssh(host=None, username=None, password=None):
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy()) 
    try:
        print "creating connection"
        ssh.connect(host, username=username, password=password)
        print "connected"
        yield ssh
    finally:
        print "closing connection"
        ssh.close()
        print "closed"    

@contextmanager
def create_sftp(ssh):
    try:
        print "creating sftp"
        sftp = ssh.open_sftp()
        print "created" 
        yield sftp
    finally:
        print "closing sftp"
        sftp.close()
        print "sftp closed"

def para(hostname,mydirectory,myfilenum,username=None,password=None,port=22):
    sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    sock.settimeout(20)
    sock.connect((hostname, port))
    my_t = paramiko.Transport(sock)
    my_t.connect(hostkey=None ,username=username, password=password, pkey=None)
    my_chan = my_t.open_session()
    my_chan.get_pty()
    my_chan.invoke_shell()
    my_sftp = paramiko.SFTP.from_transport(my_t)
    my_sftp.chdir('/users/data')
    myfilestr=os.path.join(mydirectory,myfilenum)+'.dat'
    my_sftp.get(myfilenum,myfilestr)
    my_sftp.close()
    
if __name__=="__main__":
    #
    myfilestr=r'c:\tbmno3\aug25_2009_ill\192613.dat' 
    
    if 1:
        data=reader2(myfilestr)
        pylab.imshow(data.frames[22].T)
        pylab.colorbar()
        pylab.show()
    
    
    if 0:
        myfilestr=r'c:\tbmno3\aug25_2009_ill\192547.dat' 
        data=reader(myfilestr)
        L=data.L
        y=data.intensity[:,0]
        yerr=N.sqrt(y)
        pylab.errorbar(L,y,yerr=yerr,marker='s',mfc='red',mec='red')
        pylab.semilogy()
        pylab.show()
    
    
    if 0:
        #host = '192.168.10.142'
        host='d10'
        host='192.93.249.66'
        username = 'd10'
        password = 'd10d10'
        mydirectory=r'c:\tbmno3\aug25_2009_ill'
        myfilenum='192614'
        para(host,mydirectory,myfilenum,username=username,password=password)
#    with create_ssh(host=host,username=username,password=password) as ssh:
#        with sftp_common.create_sftp(ssh) as sftp:
#            sftp.chdir('/users/data')
#            sftp.get('192553')
        
