from datetime import datetime
import time
import sys,os
import numpy as np

def get_tokenized_line(myfile,returnline=[''],splitchar=None):
    lineStr=myfile.readline()
    returnline[0]=lineStr.rstrip()
    strippedLine=lineStr.lower().rstrip().lstrip()
    if splitchar==None:
        tokenized=strippedLine.split()
    else:
        tokenized=strippedLine.split(splitchar)
    return tokenized


def fixfile(myfilestr):
    infile = open(myfilestr, "rb" )
    instr = infile.read()
    infile.close()
    outstr = instr.replace( "\r\n", "\n" ).replace( "\r", "\n" ).replace( "\n", "\r\n" )

    if len(outstr) == len(instr):
        print 'same length'
        return 

    outfile = open(myfilestr, "wb" )
    outfile.write( outstr )
    outfile.close()


def readsansfile(myfilestr):
    fixfile(myfilestr)
    myfile=open(myfilestr,"r")
    count=0
    header=[]
    myFlag=True
    data=[]
    while myFlag:
        try:
            returnline=['']
            if count <17:
                currline=get_tokenized_line(myfile,returnline)
                header.append(returnline)
                print len(currline), currline
                if count==0:
                    mytimestampstr=currline[1]+' '+currline[2]
                    mydatetime=datetime.strptime(mytimestampstr[1:-2], "%d-%b-%Y %H:%M:%S")
                    mytimestamp=time.mktime(mydatetime.timetuple())
                if count==2:
                    monitor=float(currline[1])          
            else:
                currline=get_tokenized_line(myfile,returnline,splitchar=',')
                print len(currline), currline
                if currline in ['',['']]:
                    myFlag=False
                else:
                    #data.append(np.array(currline[:-1],'Float64'))
                    data=np.concatenate((data,np.array(currline[:-1],'Float64')))
            count=count+1
        except:
            print 'done'
            myFlag=False        
    res={}
    res['header']=header
    res['monitor']=monitor
    res['timestamp']=mytimestamp
    res['data']=np.array(data)
    return res
    #Epoch=float(tokenized[1])
    #timeobj=mx.DateTime.DateTimeFromTicks(ticks=Epoch) #what I originally used
    #timeobj=datetime.fromtimestamp(Epoch)
    

def readfiles(ppfile):
    myfiledir=r"c:\bfosans\data\grasp"
    myfiledir=r"c:\bfosans\data"
    #ppfile=150
    myend='.GSP'
    myname="BFORV"+str(ppfile)
    myfilestr=os.path.join(myfiledir,myname)
    myfilestr=myfilestr+myend
    #return myfilestr
    result=readsansfile(myfilestr)
    print result
    return result
    
def read_driver(myfilestr):
    fixfile(myfilestr)
    myfile=open(myfilestr,"r")
    count=0
    myFlag=True
    res=[]
    #++    -+    +-    --    
    #* denotes absent scan
    results=[]
    while myFlag:
        try:
            currline=get_tokenized_line(myfile)
            resdict={}
            if currline[0]!='*':
                resdict['off_off']=currline[0]
            if currline[1]!='*':
                resdict['on_off']=currline[1]
            if currline[2]!='*':
                resdict['off_on']=currline[2]
            if currline[3]!='*':
                resdict['off_off']=currline[3]
            results.append(resdict)     
        except:
            print 'done'
            myFlag=False
        
    return results
            
        
        
        
    
    

if __name__=="__main__":
    myfiledir=r"c:\bfosans"
    myname="Pnki_BFORV_Scans.txt"
    myname="Ppki_BFORW_Scans.txt" # big endian??
    myname="p2.txt"
    myfilestr=os.path.join(myfiledir,myname)
    results=read_driver(myfilestr)
    keys=['off_off','off_on','on_off','on_on']
    i=0
    res={}
    for key in keys:
        try:
            fp=results[i][key]
            res[key]=readfiles(fp)
            print 'read'
        except:
            print key, 'does not exist for',i
    print res
