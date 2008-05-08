
class scanparser:
    def __init__(self,scanstr):
        self.scanstr=scanstr




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


if  __name__=='__main__':
    scanstr='Scan:Title=ICEFindPeak:Type=6:Fixed=0:FixedE=1:CountType=Time:Counts=2.0:Range=A4=50.0095 0.2:Npts=21:DetectorType=Detector:Filename=fpx:Range=A3=115.113 0.1::Title=FindPeak'