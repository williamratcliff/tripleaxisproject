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
    myfilestr=r'c:\tbmno3\aug25_2009_ill\192547.dat' 
    #host = '192.168.10.142'
    host='d10'
    host='192.93.249.66'
    username = 'd10'
    password = 'd10d10'
    mydirectory=r'c:\tbmno3\aug25_2009_ill'
#    for i in range(192533,192663):
    for i in range(192509,192511):
        #myfilenum='192601'
        myfilenum=str(i)
        para(host,mydirectory,myfilenum,username=username,password=password)
#    with create_ssh(host=host,username=username,password=password) as ssh:
#        with sftp_common.create_sftp(ssh) as sftp:
#            sftp.chdir('/users/data')
#            sftp.get('192553')
        
    sys.exit()
    data=reader(myfilestr)
    L=data.L
    y=data.intensity[:,0]
    yerr=N.sqrt(y)
    pylab.errorbar(L,y,yerr=yerr,marker='s',mfc='red',mec='red')
    pylab.semilogy()
    pylab.show()