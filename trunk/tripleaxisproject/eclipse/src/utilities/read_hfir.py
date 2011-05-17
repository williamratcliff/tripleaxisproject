import numpy as N
#import pylab
import datetime
from time import mktime
#import mx.DateTime
import writebt7
import re
import scanparser
import os
from copy import deepcopy

months={'jan':1,'feb':2,'mar':3,'apr':4,'may':5,'jun':6,'jul':7,'aug':8,'sep':9,'oct':10,'nov':11,'dec':12}

def get_tokenized_line(myfile,returnline=['']):
    lineStr=myfile.readline()
    returnline[0]=lineStr.rstrip()
    strippedLine=lineStr.lower().rstrip()
    tokenized=strippedLine.split()

    return tokenized




def get_columnmetadatas(myfile):
#get first line
    tokenized=get_tokenized_line(myfile)
    columndict={}
    columnlist=[]
    for i in N.arange(len(tokenized)):
        field=tokenized[i].lower()
        if field=='Q(x)'.lower():
            field='Qx'.lower()
        if field=='Q(y)'.lower():
            field='Qy'.lower()
        if field=='Q(z)'.lower():
            field='Qz'.lower()
        if field=='T-act'.lower():
            field='Temp'.lower()
        columndict[field]=[]
        columnlist.append(field)
    return columndict,columnlist


def readcolumns(myfile):
    columndict,columnlist=get_columnmetadatas(myfile)
    # get the names of the fields
    #prepare to read the data
    count =  0
    while 1:
        lineStr = myfile.readline()
        if not(lineStr):
            break
        if lineStr[0] != "#":
            count=count+1
            strippedLine=lineStr.rstrip().lower()
            tokenized=strippedLine.split()
            for i in range(len(tokenized)):
                field=columnlist[i]
                columndict[field].append(float(tokenized[i]))
    return columndict,columnlist

def get_columnmetadatas_bt7(tokenized):
#get first line
#    tokenized=get_tokenized_line(myfile)
    columndict={}
    columnlist=[]
    timestamp_flag=True
    #originally set the timestamp flag to True, if it turns out that there is not timestamp in the file, then create one using the time field
    for i in N.arange(1,len(tokenized)):
        field=tokenized[i]
        if field=='QX':
            field='Qx'.lower()
        if field=='QY':
            field='Qy'.lower()
        if field=='QZ':
            field='Qz'.lower()
        if field=='T-act':
            field='Temp'.lower()
        columndict[field]=[]
        columnlist.append(field)
    #In old bt7 files, there was no timestamp, so add one for those, otherwise use the one in the file
    if columndict.has_key('timestamp')==False:
        timestamp_flag=False
        columnlist.append('timestamp')
        #print 'no timestamp'
        columndict['timestamp']=[]
    return columndict,columnlist




def readfile(myfilestr,metadata):
    #get first line
    myFlag=True
    #metadata={}
    header=[]
    returnline=['']
    myfile=open(myfilestr)
    additional_metadata={}
    while myFlag:
        tokenized=get_tokenized_line(myfile,returnline=returnline)
        #print tokenized
        if tokenized==[]:
            tokenized=['']
        if tokenized[0].lower()=="#Date".lower():
            pass
        if tokenized[0].lower()=="#Date".lower():
            date_tokens=tokenized[1].split('-')
            metadata['timestamp']={}
            month=int(date_tokens[1].strip("\'"))
            day=int(date_tokens[2].strip("\'"))
            year=int(date_tokens[0].strip("\'"))
            stime=tokenized[2].strip("\'")
            stimetok=stime.split(':')
            hour=int(stimetok[0])
            minute=int(stimetok[1])
            second=int(stimetok[2])
            metadata['timestamp']['month']=int(date_tokens[1].strip("\'"))
            metadata['timestamp']['day']=int(date_tokens[2].strip("\'"))
            metadata['timestamp']['year']=int(date_tokens[0].strip("\'"))
            metadata['timestamp']['time']=tokenized[2].strip("\'")
        elif tokenized[0].lower()=="#Epoch".lower():
            #timeobj=date.datatetime(year,month,day,hour,minute,second)
            Epoch=float(tokenized[1])
            #timeobj=mx.DateTime.DateTimeFromTicks(ticks=Epoch) #what I originally used
            timeobj=datetime.datetime.fromtimestamp(Epoch)
            #print 'timeobj ',timeobj
            #print 'Epoch ', Epoch
            metadata['timestamp']['epoch']=Epoch#timeobj
            #print self.metadata['timestamp']
        elif tokenized[0].lower()=="#InstrName".lower():
            metadata['file_info']['instrument']=tokenized[1].lower()
        elif tokenized[0].lower()=="#Filename".lower():
            metadata['file_info']['filename']=tokenized[1]
            #print 'filename ', tokenized[1]
            pattern = re.compile('^(?P<base>[^.]*?)(?P<seq>[0-9]*)(?P<ext>[.].*)?$')
            match = pattern.match(tokenized[1]+'.bt7')
            dict((a,match.group(a)+"") for a in ['base','seq','ext'])
            #print 'filebase ',match.group('base')
            metadata['file_info']['filebase']=match.group('base')
            metadata['file_info']['fileseq_number']=match.group('seq')

        else:
            currfield=tokenized[0].lower().lower().strip('#')
            additional_metadata[currfield]=(tokenized[1:])
        if tokenized[1]!='Pt.'.lower():
            header.append(returnline[0])
        if tokenized[1]=='Pt.'.lower():
            columndict,columnlist=get_columnmetadatas_bt7(tokenized)
            count =  0
            try:
                lines=int(self.lines)
            except:
                lines=N.Inf
            while 1:
                lineStr = myfile.readline()
                if not(lineStr):
                    break
                if lineStr[0] != "#":
                    if count>=lines:
                        break
                    strippedLine=lineStr.rstrip()
                    tokenized=strippedLine.split()
                    for i in range(len(tokenized)):
                        field=columnlist[i]
                        columndict[field].append(float(tokenized[i]))                           
                    count=count+1
            myFlag=False
    if len(columndict[columnlist[0]])==0:
        columndict={}
        columnlist=[]
        #This is a drastic step, but if the file is empty, then no point in even recording the placeholders
    #print self.columndict['Qx']
    #print self.columnlist
    data=Data()
    data.header=deepcopy(header)
    data.data=deepcopy(columndict)
    data.metadata=deepcopy(metadata)
    data.columnlist=deepcopy(columnlist)
    data.additional_metadata=deepcopy(additional_metadata)
    return data


class Data(object):
    def __init__(self):
        self.header=[]
        self.metadata={}
        self.data={}
        self.columnlist=[]
        self.additional_metadata={}


def num2string(num):
    numstr=None
    if num<10:
        numstr='00'+str(num)
    elif (num>=10 & num <100):
        numstr='0'+str(num)
    elif (num>100):
        numstr=str(num)
    return numstr

def genfiles(myfile_nums, myfile_end='.dat',mydirectory=r'C:\hfir\HB3A\exp102\Datafiles',myfile_base='HB3A_exp0102_scan'):
    file_list=[]
    for myfile_num in myfile_nums:
        if myfile_num<10:
            myfile_numstr='000'+str(myfile_num)
        elif myfile_num>9 and myfile_num<100:
            myfile_numstr='00'+str(myfile_num)
        elif myfile_num>99 and myfile_num<1000:
            myfile_numstr='0'+str(myfile_num)
        else:
            myfile_numstr=str(myfile_num)
        myfilestr=os.path.join(mydirectory,myfile_base+myfile_numstr+myfile_end)
        print myfilestr
        file_list.append(myfilestr)
    return file_list

if __name__=='__main__':
    if 1:
        mon0=9000
        myfilestr=r'c:\bifeo3xtal\jan8_2008\9175\mesh53439.bt7'
       
        
        
        #0014.dat
        myfile_num=84
        file_list=genfiles(range(84,100))
        print file_list
        #metadata={}
        #data=readfile(myfilestr,metadata)
        #print data.header
        #print data.metadata
        #print data.columnlist



