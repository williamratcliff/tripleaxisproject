import numpy as N
class scanparser:
    def __init__(self,scanstr):
        self.scanstr=scanstr
        self.scan_description={}
    def parse_range(self,rangestr):
        print rangestr
        prange={}
        npts=self.scan_description['npts']
        print 'npts',npts
        range_split=rangestr.split('range=')
        print 'range_split',range_split
        fields=range_split[1].split('=')
        #print fields
        if fields[0]=='e':
            prange['e']={}
            toks=fields[1].split()
            if len(toks)==3:
                if toks[-1]=='s':
                    print 'start stop'
                    prange['e']['start']=min(float(toks[0]),float(toks[1]))
                    prange['e']['stop']=max(float(toks[0]),float(toks[1]))
                    if npts>1:
                        prange['e']['step']=float(prange['e']['stop']-prange['e']['start'])/(npts-1)
                    else:
                        prange['e']['step']=float(0)
                elif toks[-1]=='i':
                    print 'start increment'
                    start=float(toks[0])
                    step=float(toks[1])
                    stop=start+(npts-1)*step
                    prange['e']['step']=N.absolute(step)
                    prange['e']['start']=min(start,stop)
                    prange['e']['stop']=max(start,stop)

            else:
                print 'center step'
                step=float(toks[1])
                center=float(toks[0])
                print 'center',center
                start=center-float(step)*(npts-1)/2
                stop=center+float(step)*(npts-1)/2
                prange['e']['step']=N.absolute(step)
                prange['e']['start']=min(start,stop)
                prange['e']['stop']=max(start,stop)
        if 1:
            if fields[0]=='q':
                toks=fields[1].split()
                prange['q']={}
                #print 'toks',toks
                if len(toks)==3:
                    if toks[-1]=='s':
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
        return prange


    def parse_scan(self):
        scanstr=self.scanstr
        self.scan_description={}
        scan_description=self.scan_description
        scan_description['scan_string']=scanstr
        scan_description['range_strings']=[]

        toks=scanstr.split(':')
        try:
            if toks[0].lower()!='scan':
                raise BadScanError,'Not a Valid Scan'
            toks=toks[1:]
            for tok in toks:
                field=tok.split('=')
                key=field[0].lower()
                value=field[1]
                if key.lower()=='range':
                    scan_description['range_strings'].append(tok.lower())
                else:
                    try:
                        scan_description[key]=float(value)
                    except ValueError:
                        scan_description[key]=value
            return self.scan_description
        except BadScanError:
            print 'Not a Valid Scan'
            return None

class BadScanError(Exception):
     def __init__(self, value):
         self.value = value
     def __str__(self):
         return repr(self.value)

if  __name__=='__main__':
#find peak, A3-A4
    if 1:
        scanstr='Scan:Title=ICEFindPeak:Type=6:Fixed=0:FixedE=1:CountType=Time:Counts=2.0:Range=A4=50.0095 0.2:Npts=21:DetectorType=Detector:Filename=fpx:Range=A3=115.113 0.1::Title=FindPeak'
#inititial final h
    if 1:
        scanstr='Scan:SubID=13176:JType=VECTOR:Fixed=1:FixedE=13.6998911684:Npts=1:Counts=1.0:Prefac=1.0:DetectorType=Detector:CountType=Time:Filename=dumb:HoldScan=0.0:Range=Q=1.0~0.0~0.0 2.0~0.0~0.0 s:Range=E=0.0 0.0 s'
#initial step h
    if 1:
        scanstr='Scan:SubID=13176:JType=VECTOR:Fixed=1:FixedE=13.6998911684:Npts=1:Counts=1.0:Prefac=1.0:DetectorType=Detector:CountType=Time:Filename=dumb:HoldScan=0.0:Range=Q=1.0~0.0~0.0 2.0~0.0~0.0 i:Range=E=0.0 0.0 i'
#center step h
    if 1:
        scanstr='Scan:SubID=13176:JType=VECTOR:Fixed=1:FixedE=13.6998911684:Npts=1:Counts=1.0:Prefac=1.0:DetectorType=Detector:CountType=Time:Filename=dumb:HoldScan=0.0:Range=Q=1.0~0.0~0.0 2.0~0.0~0.0:Range=E=0.0 0.0'

#center step e [-1,0,1]
    if 0:
        scanstr='Scan:SubID=13176:JType=VECTOR:Fixed=1:FixedE=13.6998911684:Npts=3:Counts=1.0:Prefac=1.0:DetectorType=Detector:CountType=Monitor:Filename=dumb:HoldScan=0.0:Range=Q=0.0~0.0~0.0 0.0~0.0~0.0:Range=E=0.0 1.0'
#center step e [-.5,.5]
    if 0:
        scanstr='Scan:SubID=13176:JType=VECTOR:Fixed=1:FixedE=13.6998911684:Npts=2:Counts=1.0:Prefac=1.0:DetectorType=Detector:CountType=Monitor:Filename=dumb:HoldScan=0.0:Range=Q=0.0~0.0~0.0 0.0~0.0~0.0:Range=E=0.0 1.0'
#initial step e [0,1,2]
    if 0:
        scanstr='Scan:SubID=13176:JType=VECTOR:Fixed=1:FixedE=13.6998911684:Npts=3:Counts=1.0:Prefac=1.0:DetectorType=Detector:CountType=Monitor:Filename=dumb:HoldScan=0.0:Range=Q=0.0~0.0~0.0 0.0~0.0~0.0 i:Range=E=0.0 1.0 i'
#start stop e [0,.5,1]
    if 0:
        scanstr='Scan:SubID=13176:JType=VECTOR:Fixed=1:FixedE=13.6998911684:Npts=3:Counts=1.0:Prefac=1.0:DetectorType=Detector:CountType=Monitor:Filename=dumb:HoldScan=0.0:Range=Q=0.0~0.0~0.0 0.0~0.0~0.0 s:Range=E=0.0 1.0 s'
    myparser=scanparser(scanstr)
    scanstr_parsed=myparser.parse_scan()
    print myparser.parse_range(scanstr_parsed['range_strings'][1])