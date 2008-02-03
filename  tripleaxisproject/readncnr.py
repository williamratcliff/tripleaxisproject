import numpy as N
import  pylab
import datetime
import mx.DateTime
import writebt7


months={'jan':1,'feb':2,'mar':3,'apr':4,'may':5,'jun':6,'jul':7,'aug':8,'sep':9,'oct':10,'nov':11,'dec':12}

def get_tokenized_line(myfile,returnline=['']):
        lineStr=myfile.readline()
        returnline[0]=lineStr.rstrip()
        strippedLine=lineStr.lower().rstrip()
        tokenized=strippedLine.split()
        
        return tokenized


class datareader:
    def __init__(self,myfilestr=None):
        self.myfilestr=myfilestr
        #define Data Abstraction Layer
        self.metadata={}
        self.metadata['count_info']={}
        self.metadata['count_info']['monitor_base']=None #float(tokenized[6])
        self.metadata['count_info']['monitor_prefactor']=None#float(tokenized[7])
        self.metadata['count_info']['monitor']=None#self.metadata['count_info']['monitor_base']*self.metadata['count_info']['monitor_prefactor']
        self.metadata['count_info']['count_type']=None  #can be 'monitor', 'time' #tokenized[8].strip("'").lower()
                
        self.metadata['file_info']={}
        self.metadata['file_info']['filename']=None#tokenized[0].strip("'")
        self.metadata['file_info']['filebase']=None#self.metadata['file_info']['filename'][0:5]
        self.metadata['file_info']['scantype']=None#tokenized[5].strip("'").lower()
        self.metadata['file_info']['instrument']=None#self.metadata['file_info']['filename'].split('.')[1].lower()
        self.metadata['file_info']['comment']=None #myfile.readline().rstrip()
        self.metadata['file_info']['scan_description']=None
        self.metadata['file_info']['experiment_id']=None
        
        
        self.metadata['timestamp']={}
        self.metadata['timestamp']['month']=None#int, for icp data it is translated using the months dict
        self.metadata['timestamp']['day']=None#int
        self.metadata['timestamp']['year']=None#int
        self.metadata['timestamp']['time']=None#str

        self.metadata['collimations']={}
        self.metadata['collimations']['coll1']=None#float(tokenized[0])
        self.metadata['collimations']['coll2']=None#float(tokenized[1])
        self.metadata['collimations']['coll3']=None#float(tokenized[2])
        self.metadata['collimations']['coll4']=None#float(tokenized[3])

        self.metadata['mosaic']={}
        self.metadata['mosaic']['mosaic_monochromator']=None#float(tokenized[4])
        self.metadata['mosaic']['mosaic_sample']=None#float(tokenized[5])
        self.metadata['mosaic']['mosaic_analyzer']=None#float(tokenized[6])
        
        self.metadata['energy_info']={}
        self.metadata['energy_info']['wavelength']=float(tokenized[7])
        
        self.metadata['temperature_info']={}
        self.metadata['temperature_info']['Tstart']=float(tokenized[8])
        self.metadata['temperature_info']['Tstep']=float(tokenized[9])
        
        self.metadata['magnetic_field']={}
        self.metadata['magnetic_field']['Hfield']=float(tokenized[10])

        self.metadata['orient1']={}
        self.metadata['orient1']['h']=float(tokenized[7])
        self.metadata['orient1']['k']=float(tokenized[8])
        self.metadata['orient1']['l']=float(tokenized[9])
        #ignore "angle" field
        self.metadata['orient2']={}
        self.metadata['orient2']['h']=float(tokenized[11])
        self.metadata['orient2']['k']=float(tokenized[12])
        self.metadata['orient2']['l']=float(tokenized[13])

        #self.metadata['lattice']={}
        self.metadata['lattice']['a']=float(tokenized[0])
        self.metadata['lattice']['b']=float(tokenized[1])
        self.metadata['lattice']['c']=float(tokenized[2])
        self.metadata['lattice']['alpha']=float(tokenized[3])
        self.metadata['lattice']['beta']=float(tokenized[4])
        self.metadata['lattice']['gamma']=float(tokenized[5])

        self.metadata['energy_info']={}
        self.metadata['q_center']={}
        self.metadata['q_step']={}
        self.metadata['q_center']['e_center']=float(tokenized[0])
        self.metadata['q_step']['delta_e']=float(tokenized[1])
        self.metadata['energy_info']['ef']=float(tokenized[2])
        
        self.metadata['dspacing']={}
        self.metadata['dspacing']['monochromator_dspacing']=float(tokenized[3])
        self.metadata['dspacing']['analyzer_dspacing']=float(tokenized[4])
        
        self.metadata['temperature_info']={}
        self.metadata['temperature_info']['Tstart']=float(tokenized[5])
        self.metadata['temperature_info']['Tstep']=float(tokenized[6])
        tokenized=get_tokenized_line(myfile)
        self.metadata['energy_info']['efixed']=tokenized[4]
       
        self.metadata['q_center']['h_center']=float(tokenized[0])
        self.metadata['q_center']['k_center']=float(tokenized[1])
        self.metadata['q_center']['l_center']=float(tokenized[2])
        self.metadata['q_step']['delta_h']=float(tokenized[3])
        self.metadata['q_step']['delta_k']=float(tokenized[4])
        self.metadata['q_step']['delta_l']=float(tokenized[5])
        self.metadata['magnetic_field']['hfield']=float(tokenized[6])


    def readimotors(self,myfile):
    #motor1
        tokenized=get_tokenized_line(myfile)
    #    print tokenized
        motor1={'start':float(tokenized[1])}
        motor1['step']=float(tokenized[2])
        motor1['end']=float(tokenized[3])
        self.metadata['motor1']=motor1

    #motor2
        tokenized=get_tokenized_line(myfile)
    #    print tokenized
        motor2={'start':float(tokenized[1])}
        motor2['step']=float(tokenized[2])
        motor2['end']=float(tokenized[3])
        self.metadata['motor2']=motor2
        
    #motor3
        tokenized=get_tokenized_line(myfile)
    #    print tokenized
        motor3={'start':float(tokenized[1])}
        motor3['step']=float(tokenized[2])
        motor3['end']=float(tokenized[3])
        self.metadata['motor3']=motor3

    #motor4
        tokenized=get_tokenized_line(myfile)
    #    print tokenized
        motor4={'start':float(tokenized[1])}
        motor4['step']=float(tokenized[2])
        motor4['end']=float(tokenized[3])
        self.metadata['motor4']=motor4

    #motor5
        tokenized=get_tokenized_line(myfile)
    #    print tokenized
        motor5={'start':float(tokenized[1])}
        motor5['step']=float(tokenized[2])
        motor5['end']=float(tokenized[3])
        self.metadata['motor5']=motor5

    #motor6
        tokenized=get_tokenized_line(myfile)
    #    print tokenized
        motor6={'start':float(tokenized[1])}
        motor6['step']=float(tokenized[2])
        motor6['end']=float(tokenized[3])
        self.metadata['motor6']=motor6
        #skip line describing Motor Start Step End
        lineStr = myfile.readline()
        return
    
    def readimetadata(self,myfile):
    #experiment info
        tokenized=get_tokenized_line(myfile)
        #collimations=[] #in stream order
        
        #self.metadata['collimations']={}
        self.metadata['collimations']['coll1']=float(tokenized[0])
        self.metadata['collimations']['coll2']=float(tokenized[1])
        self.metadata['collimations']['coll3']=float(tokenized[2])
        self.metadata['collimations']['coll4']=float(tokenized[3])
        #collimations.append(float(tokenized[1]))
        #collimations.append(float(tokenized[2]))
        #collimations.append(float(tokenized[3]))
        
        #mosaic=[] #order is monochromator, sample, mosaic
        #self.metadata['mosaic']={}
        self.metadata['mosaic']['mosaic_monochromator']=float(tokenized[4])
        self.metadata['mosaic']['mosaic_sample']=float(tokenized[5])
        self.metadata['mosaic']['mosaic_analyzer']=float(tokenized[6])
        
        #self.metadata['energy_info']={}
        self.metadata['energy_info']['wavelength']=float(tokenized[7])
        
        #self.metadata['temperature_info']={}
        self.metadata['temperature_info']['Tstart']=float(tokenized[8])
        self.metadata['temperature_info']['Tstep']=float(tokenized[9])
        
        #self.metadata['magnetic_field']={}
        self.metadata['magnetic_field']['Hfield']=float(tokenized[10])
        #print tokenized
        #skip field names of experiment info
        lineStr=myfile.readline()
        self.readimotors(myfile)
        return


    def readqmetadata(self,myfile):
        #experiment info
        tokenized=get_tokenized_line(myfile)
##        collimations=[] #in stream order
##        collimations.append(float(tokenized[0]))
##        collimations.append(float(tokenized[1]))
##        collimations.append(float(tokenized[2]))
##        collimations.append(float(tokenized[3]))
##        self.metadata['collimations']=collimations
##        mosaic=[] #order is monochromator, sample, mosaic
##        mosaic.append(float(tokenized[4]))
##        mosaic.append(float(tokenized[5]))
##        mosaic.append(float(tokenized[6]))
##        self.metadata['mosaic']=mosaic

        
        #self.metadata['collimations']={}
        self.metadata['collimations']['coll1']=float(tokenized[0])
        self.metadata['collimations']['coll2']=float(tokenized[1])
        self.metadata['collimations']['coll3']=float(tokenized[2])
        self.metadata['collimations']['coll4']=float(tokenized[3])
        
        #self.metadata['mosaic']={}
        self.metadata['mosaic']['mosaic_monochromator']=float(tokenized[4])
        self.metadata['mosaic']['mosaic_sample']=float(tokenized[5])
        self.metadata['mosaic']['mosaic_analyzer']=float(tokenized[6])

        
        #self.metadata['orient1']={}
        self.metadata['orient1']['h']=float(tokenized[7])
        self.metadata['orient1']['k']=float(tokenized[8])
        self.metadata['orient1']['l']=float(tokenized[9])
        #ignore "angle" field
        #self.metadata['orient2']={}
        self.metadata['orient2']['h']=float(tokenized[11])
        self.metadata['orient2']['k']=float(tokenized[12])
        self.metadata['orient2']['l']=float(tokenized[13])
        
##        orient1.append(float(tokenized[7]))
##        orient1.append(float(tokenized[8]))
##        orient1.append(float(tokenized[9]))
##        self.metadata['orient1']=orient1
##        #ignore the "angle" field
##        orient2=[]
##        orient2.append(float(tokenized[11]))
##        orient2.append(float(tokenized[12]))
##        orient2.append(float(tokenized[13]))
##        self.metadata['orient2']=orient2
        #skip line with field names
        myfile.readline()
        tokenized=get_tokenized_line(myfile)
        #self.metadata['lattice']={}
        self.metadata['lattice']['a']=float(tokenized[0])
        self.metadata['lattice']['b']=float(tokenized[1])
        self.metadata['lattice']['c']=float(tokenized[2])
        self.metadata['lattice']['alpha']=float(tokenized[3])
        self.metadata['lattice']['beta']=float(tokenized[4])
        self.metadata['lattice']['gamma']=float(tokenized[5])
        #self.metadata['lattice']=lattice
        #skip line with field names
        myfile.readline()
        tokenized=get_tokenized_line(myfile)
        #self.metadata['energy_info']={}
        #self.metadata['q_center']={}
        #self.metadata['q_step']={}
        self.metadata['q_center']['e_center']=float(tokenized[0])
        self.metadata['q_step']['delta_e']=float(tokenized[1])
        self.metadata['energy_info']['ef']=float(tokenized[2])
        #self.metadata['dspacing']={}
        self.metadata['dspacing']['monochromator_dspacing']=float(tokenized[3])
        self.metadata['dspacing']['analyzer_dspacing']=float(tokenized[4])
        #self.metadata['temperature_info']={}
        self.metadata['temperature_info']['Tstart']=float(tokenized[5])
        self.metadata['temperature_info']['Tstep']=float(tokenized[6])
        tokenized=get_tokenized_line(myfile)
        self.metadata['energy_info']['efixed']=tokenized[4]
        tokenized=get_tokenized_line(myfile)
        #qcenter=[]
        #qstep=[]
        self.metadata['q_center']['h_center']=float(tokenized[0])
        self.metadata['q_center']['k_center']=float(tokenized[1])
        self.metadata['q_center']['l_center']=float(tokenized[2])
        self.metadata['q_step']['delta_h']=float(tokenized[3])
        self.metadata['q_step']['delta_k']=float(tokenized[4])
        self.metadata['q_step']['delta_l']=float(tokenized[5])
        #self.metadata['magnetic_field']={}
        self.metadata['magnetic_field']['hfield']=float(tokenized[6])
        #skip line describing fields
        linestr=myfile.readline()
        return

    def readbmetadata(self,myfile):
        self.readqmetadata(myfile)
        self.readimotors(myfile)
        return



    def get_columnmetadatas(self,myfile):
    #get first line
        tokenized=get_tokenized_line(myfile)
        self.columndict={}
        #self.columndict['columnlist']=[]
        self.columnlist=[]
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
            self.columndict[field]=[]
            #self.columndict['columnlist'].append(field)
            self.columnlist.append(field)
        return 

    def determinefiletype(self,myfile):
    #get first line
        tokenized=get_tokenized_line(myfile)
        #self.metadata={'filetype':tokenized[5].strip("'")}
        #self.metadata={}
        #self.metadata['count_info']={}
        self.metadata['count_info']['monitor_base']=float(tokenized[6])
        self.metadata['count_info']['monitor_prefactor']=float(tokenized[7])
        self.metadata['count_info']['monitor']=self.metadata['count_info']['monitor_base']*self.metadata['count_info']['monitor_prefactor']
        self.metadata['count_info']['count_type']=tokenized[8].strip("'").lower()
                
        #self.metadata['file_info']={}
        self.metadata['file_info']['filename']=tokenized[0].strip("'")
        self.metadata['file_info']['filebase']=self.metadata['file_info']['filename'][0:5]
        self.metadata['file_info']['scantype']=tokenized[5].strip("'").lower()
        self.metadata['file_info']['instrument']=self.metadata['file_info']['filename'].split('.')[1].lower()
        
        #self.metadata['timestamp']={}
        month_str=tokenized[1].strip("\'").lower()
        
        self.metadata['timestamp']['month']=months[month_str]#tokenized[1].strip("\'").lower()
        self.metadata['timestamp']['day']=int(tokenized[2].strip("\'"))
        self.metadata['timestamp']['year']=int(tokenized[3].strip("\'"))
        self.metadata['timestamp']['time']=tokenized[4].strip("\'")
        
        #I put this away for now, because it is not reliable about the actual number of points in the file, just the desired number
        #self.metadata['npts']=int(tokenized[9])
        
        
        
        #skip over names of fields 
        lineStr=myfile.readline()
        #comment and filename
        self.metadata['file_info']['comment']=myfile.readline().rstrip()
        return self.metadata['file_info']['scantype']

    def readcolumns(self,myfile):
        self.get_columnmetadatas(myfile)
        # get the names of the fields
    #   prepare to read the data    
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
                    field=self.columnlist[i]
                    self.columndict[field].append(float(tokenized[i]))
        return

    def get_columnmetadatas_bt7(self,tokenized):
    #get first line
    #    tokenized=get_tokenized_line(myfile)
        self.columndict={}
        #self.columndict['columnlist']=[]
        self.columnlist=[]
        for i in N.arange(1,len(tokenized)):
            field=tokenized[i]
            if field=='QX':
                field='Qx'
            if field=='QY':
                field='Qy'
            if field=='QZ':
                field='Qz'
            if field=='T-act':
                field='Temp'
            self.columndict[field]=[]
            #self.columndict['columnlist'].append(field)
            self.columnlist.append(field)
        self.columnlist.append('timestamp')
        self.columndict['timestamp']=[]
        return 
    
    def range_parser(self,rangestr,scan_description):
        #print rangestr
        range_split=rangestr.split('range=')
        fields=range_split[1].split('=')
        #print fields
        if fields[0]=='e':
            toks=fields[1].split()
            if len(toks)==3:
                scan_description['range']['e']['start']=float(toks[0])
                scan_description['range']['e']['stop']=float(toks[1])
            else:
                scan_description['range']['e']['center']=float(toks[0])
                scan_description['range']['e']['step']=float(toks[1])
        
        if fields[0]=='q':
            toks=fields[1].split()
            #print 'toks',toks
            if len(toks)==3:
                start=toks[0].split('~')
                scan_description['range']['q']['start']={}
                scan_description['range']['q']['stop']={}
                scan_description['range']['q']['start']['h']=float(start[0])
                scan_description['range']['q']['start']['k']=float(start[1])
                scan_description['range']['q']['start']['l']=float(start[2])
                stop=toks[0].split('~')
                scan_description['range']['q']['stop']['h']=float(stop[0])
                scan_description['range']['q']['stop']['k']=float(stop[1])
                scan_description['range']['q']['stop']['l']=float(stop[2])
            else:
                start=toks[0].split('~')
                scan_description['range']['q']['center']={}
                scan_description['range']['q']['step']={}
                scan_description['range']['q']['center']['h']=float(start[0])
                scan_description['range']['q']['center']['k']=float(start[1])
                scan_description['range']['q']['center']['l']=float(start[2])
                stop=toks[0].split('~')
                scan_description['range']['q']['step']['h']=float(stop[0])
                scan_description['range']['q']['step']['k']=float(stop[1])
                scan_description['range']['q']['step']['l']=float(stop[2])
        
        #print 'Range', scan_description['range']
        return
        
    
    def parse_scan(self,scanstr):
        scan_description={}
        scan_description['scan_string']=scanstr
        scan_description['range']={}
        scan_description['range']['e']={}
        scan_description['range']['q']={}
        
        toks=scanstr.split(':')
        for i in range(1,len(toks)):
            if toks[i][0]=='r':
                self.range_parser(toks[i],scan_description)
            else:
                fields=toks[i].split('=')
                try:
                    scan_description[fields[0]]=float(fields[1])
                except ValueError:
                    scan_description[fields[0]]=(fields[1])
        return scan_description

    def readbt7(self,myfile):
    #get first line
        myFlag=True
        self.metadata={}
        self.header=[]
        self.metadata['dspacing']={}
        returnline=['']
        while myFlag:
            tokenized=get_tokenized_line(myfile,returnline=returnline)
            #print tokenized
            if tokenized[0].lower()=="#Date".lower():
                date_tokens=tokenized[1].split('-')
                self.metadata['timestamp']={}
                month=int(date_tokens[1].strip("\'"))
                day=int(date_tokens[2].strip("\'"))
                year=int(date_tokens[0].strip("\'"))
                stime=tokenized[2].strip("\'")
                stimetok=stime.split(':')
                hour=int(stimetok[0])
                minute=int(stimetok[1])
                second=int(stimetok[2])
                self.metadata['timestamp']['month']=int(date_tokens[1].strip("\'"))
                self.metadata['timestamp']['day']=int(date_tokens[2].strip("\'"))
                self.metadata['timestamp']['year']=int(date_tokens[0].strip("\'"))
                self.metadata['timestamp']['time']=tokenized[2].strip("\'")
            elif tokenized[0].lower()=="#Epoch".lower(): 
                #timeobj=date.datatetime(year,month,day,hour,minute,second)
                Epoch=float(tokenized[1])
                timeobj=mx.DateTime.DateTimeFromTicks(ticks=Epoch)
                self.metadata['timestamp']['Epoch']=timeobj
                #print self.metadata['timestamp']
            elif tokenized[0].lower()=="#MonoSpacing".lower(): 
                #self.metadata['dmono']=float(tokenized[1])
                self.metadata['dspacing']['monochromator_dspacing']=float(tokenized[1])
            elif tokenized[0].lower()=="#AnaSpacing".lower():
                #self.metadata['dana']=float(tokenized[1]) 
                self.metadata['dspacing']['analyzer_dspacing']=float(tokenized[1])
            elif tokenized[0].lower()=="#Orient".lower():
                #self.metadata['orient1']={}
                self.metadata['orient1']['h']=float(tokenized[1])
                self.metadata['orient1']['k']=float(tokenized[2])
                self.metadata['orient1']['l']=float(tokenized[3])
                #self.metadata['orient2']={}
                self.metadata['orient2']['h']=float(tokenized[4])
                self.metadata['orient2']['k']=float(tokenized[5])
                self.metadata['orient2']['l']=float(tokenized[6])
            elif tokenized[0].lower()=="#Lattice".lower():
                #self.metadata['lattice']={}
                self.metadata['lattice']['a']=float(tokenized[1])
                self.metadata['lattice']['b']=float(tokenized[2])
                self.metadata['lattice']['c']=float(tokenized[3])
                self.metadata['lattice']['alpha']=float(tokenized[4])
                self.metadata['lattice']['beta']=float(tokenized[5])
                self.metadata['lattice']['gamma']=float(tokenized[6])
            elif tokenized[0].lower()=="#AnalyzerDetectorMode".lower():
                self.metadata['analyzerdetectormode']=tokenized[2].lower()
            elif tokenized[0].lower()=="#Signal".lower():
                self.metadata['signal']=tokenized[2].lower()
                #print self.metadata['signal']
            elif tokenized[0].lower()=="#AnalyzerDetectorDevicesOfInterest".lower():
                self.metadata['#AnalyzerDetectorDevicesOfInterest'.lower()]=tokenized[1:]
                #print self.metadata['#AnalyzerDetectorDevicesOfInterest'.lower()]
            elif tokenized[0].lower()=="#ScanDescr".lower():
                #self.metadata['scan_description']={}
                scanstr=''
                for i in range(1,len(tokenized)):
                    scanstr=scanstr+tokenized[i]+' '
                #print 'reached'
                self.metadata['scan_description']=self.parse_scan(scanstr)
                #print self.metadata['scan_description']['range']
            else:
                currfield=tokenized[0].lower().lower().strip('#')
                self.metadata[currfield]=(tokenized[1:])    
            if tokenized[0]!='#Columns'.lower():
                self.header.append(returnline[0])
            if tokenized[0]=='#Columns'.lower():
                self.get_columnmetadatas_bt7(tokenized)
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
                            field=self.columnlist[i]
                            try:
                                if field.lower()=='time':
                                    timedelta=mx.DateTime.DateTimeDelta(0,0,0,float(tokenized[i]))
                                    self.columndict['timestamp'].append((timeobj+timedelta).ticks())
                                self.columndict[field].append(float(tokenized[i]))
                            except ValueError:
                                self.columndict[field].append((tokenized[i]))
                        count=count+1
                myFlag=False
        #print self.columndict['Qx']
        #print self.columnlist
        return



    def readbuffer(self,myfilestr,lines=N.Inf):
        self.myfilestr=myfilestr
        self.lines=lines
        myfile = open(myfilestr, 'r')
        self.instrument=myfilestr.split('.')[1]
        if self.instrument in ['bt9','ng5']:
            # Determine FileType
            self.determinefiletype(myfile)
            if self.metadata['file_info']['scantype'].lower()=='i':
                print "calling readibuffer"
                self.readimetadata(myfile)
            if self.metadata['file_info']['scantype'].lower()=='b':
                print "calling readbbuffer"
                self.readbmetadata(myfile)
            if self.metadata['file_info']['scantype'].lower()=='q':
                print "calling readqbuffer"
                self.readqmetadata(myfile)
            
            #read columns
            self.readcolumns(myfile)
            myfile.close()
            mydata=Data(self.metadata,self.columndict)
            #print self.metadata
            #print self.columnlist   
        else:
            #instrument is bt7
            self.readbt7(myfile)
            #self.readbt7columns(myfile)
            myfile.close()
            if self.header==None:
                self.header=[]
            mydata=Data(self.metadata,self.columndict,self.header)
        return mydata


class Data:
    def __init__(self,metadata,data,header=None):
        self.metadata=metadata
        self.data=data
        self.header=header
    
    def get_monitor(self):
        return self.metadata['monitor']
    #@property
    #def monitor(self):
    #    "The monitor rate"
    #    def fget(self):
    #        return self.metadata['monitor']
    ##def get_filetype(self):
    ##    return self.metadata['filetype']
    #def get_data_fields(self):
    #    return self.data['columnlist']
    def get_motor1(self):
        return self.metadata['motor1']
    def get_motor2(self):
        return self.metadata['motor2']
    def get_motor3(self):
        return self.metadata['motor3']
    def get_motor4(self):
        return self.metadata['motor4']
    def get_motor5(self):
        return self.metadata['motor5']
    def get_motor6(self):
        return self.metadata['motor6']
    def get_field(self,field):
        return self.data[field]
    def gen_motor1_arr(self):
        motor=self.get_motor1()
        step=motor['step']
        start=motor['start']
        if step==0.0:
            res=start*N.ones((1,self.npts),'d')
        else:
            res=N.arange(start,motor['end'],step)
        return res
    def gen_motor2_arr(self):
        motor=self.get_motor2()
        step=motor['step']
        start=motor['start']
        if step==0.0:
            res=start*N.ones((1,self.npts),'d')
        else:
            res=N.arange(start,motor['end'],step)
        return res
    def gen_motor3_arr(self):
        motor=self.get_motor3()
        step=motor['step']
        start=motor['start']
        if step==0.0:
            res=start*N.ones((1,self.npts),'d')
        else:
            res=N.arange(start,motor['end'],step)
        return res
    def gen_motor4_arr(self):
        motor=self.get_motor4()
        step=motor['step']
        start=motor['start']
        if step==0.0:
            res=start*N.ones((1,self.npts),'d')
        else:
            res=N.arange(start,motor['end'],step)
        return res
    def gen_motor5_arr(self):
        motor=self.get_motor5()
        step=motor['step']
        start=motor['start']
        if step==0.0:
            res=start*N.ones((1,self.npts),'d')
        else:
            res=N.arange(start,motor['end'],step)
        return res        
    def gen_motor6_arr(self):
        motor=self.get_motor6()
        step=motor['step']
        start=motor['start']
        if step==0.0:
            res=start*N.ones((1,self.npts),'d')
        else:
            res=N.arange(start,motor['end'],step)
        return res
    

#   self.columndict[field]
    
    #count_type=property(get_count_type)
    #filetype=property(get_filetype)
    #npts=property(get_npts)
    motor1=property(get_motor1)
    motor2=property(get_motor2)
    motor3=property(get_motor3)
    motor4=property(get_motor4)
    motor5=property(get_motor5)
    motor6=property(get_motor6)    
    #data_fields=property(get_data_fields)
    #monitor=property(get_monitor)
    
class DataCollection:
    def __init__(self):
        self.data=[]
    def get_data(self):
        return self.data
    def add_datum(self,datum):
        self.data.append(datum)
        return
    def extract_a4(self):
        a4=[]
        for i in range(len(self.data)):
            motor4=self.data[i].get_motor4()
            a4.append(motor4['start'])
        return N.array(a4,'d')
    def extract_a3a4(self):
        a4=[]
        a3=[]
        counts=[]
        for i in range(len(self.data)):
            motor4=self.data[i].get_motor4()
            a4.append(self.data[i].gen_motor4_arr())
            a3.append(self.data[i].get_field('A3'))
            counts.append(self.data[i].get_field('COUNTS'))
        return N.ravel(N.array(a3,'d')),N.ravel(N.array(a4,'d')),N.ravel(N.array(counts,'d'))
    
    data=property(get_data,add_datum)

def num2string(num):
    numstr=None
    if num<10:
        numstr='00'+str(num)
    elif (num>=10 & num <100):
        numstr='0'+str(num)
    elif (num>100):
        numstr=str(num)
    return numstr

if __name__=='__main__':

    if 0:
        #ibuff
        myfilestr=r'c:\summerschool2007\\qCdCr014.ng5'
    if 1:
        myfilestr=r'c:\bifeo3xtal\jan8_2008\9175\meshbefieldneg1p3plusminus53470.bt7'
        #myfilestr=r'c:\bifeo3xtal\jan8_2008\9175\mesh53439.bt7'
        mydatareader=datareader()
        mydata=mydatareader.readbuffer(myfilestr,lines=91)
        myoutfilestr=r'c:\bifeo3xtal\jan8_2008\9175\meshbefieldneg1p3plusminus53470.bt7.out'
        mywriter=writebt7.datawriter()
        mywriter.write(myoutfilestr=myoutfilestr,mydata=mydata) 
        print 'done'          
        mydataout=mydata=mydatareader.readbuffer(myoutfilestr,lines=91)
        print N.array(mydata.data['qy'])-N.array(mydataout.data['qy']) 
        #print len(mydata.data['timestamp'])
        #print mydata.data['Qy']
        #print mydata.data
        #print mydata.metadata
    if 0:
        #bragg
        myfilestr=r'c:\sqltest\\nuc10014.bt9'
    if 0:
        #qbuff
        myfilestr=r'c:\sqltest\\mnl1p004.ng5'
    if 0:
        mydirectory=r'c:\summerschool2007\\' 
        myfilenumbers=range(4,33,1)
        myend='.ng5'
        myfilehead='qCdCr'
        #myfilestr=mydirectory+'qCdCr006.ng5'
        data=DataCollection()
        mydatareader=datareader()
        mydata=mydatareader.readbuffer(myfilestr)
        print mydata.metadata
#    print mydata.npts
#    print mydata.monitor
#    print mydata.gen_motor6_arr()
#    print mydata.gen_motor5_arr()
#    print mydata.gen_motor4_arr()
#    print mydata.gen_motor3_arr()
#    print mydata.gen_motor2_arr()
#    print mydata.gen_motor1_arr()
    if 0:
        for i in range(len(myfilenumbers)):
            myfilenum=num2string(myfilenumbers[i])
            myfilestr=mydirectory+myfilehead+myfilenum+myend
            print myfilestr
            data.add_datum(mydatareader.readibuffer(myfilestr))
        a3,a4,counts=data.extract_a3a4()
        print a3.shape
        print a4.shape
        print counts.shape




