class datawriter:
    def __init__(self,mydata=None,myoutfilestr=None):
        self.mydata=mydata
        self.myoutfilestr=myoutfilestr
        return
    
    def write(self,mydata=None,myoutfilestr=None):
#        myoutfilestr=r'c:\bifeo3xtal\jan8_2008\9175\meshbefieldneg1p3plusminus53470.bt7.out'      
        if mydata==None:
            mydata=self.mydata
        if myoutfilestr==None:
            myoutfilestr=self.myoutfilestr  
        count=1
        #mydata.additional_metadata['parsed_scandescription'] #TODO choose correct field based on this
        for key in mydata.data.keys():
            if key=='detector':
                detectorpos=count
                print 'detectorpos ',detectorpos
            if key=='qx':
                scanpos=count
            count=count+1
        myoutfile=open(myoutfilestr,'wt')
        for i in range(len(mydata.header)):
            s=mydata.header[i]+'\n'
            tokenized=s.rstrip().lower().split()
            if tokenized[0]=='#signal'.lower():
                s='#signal'+' '+str(detectorpos)+' '+'detector\n'
            if tokenized[0]=='#scan'.lower():
                s='#scan'+' '+str(scanpos)+' '+'qx\n'


                    
                    
            myoutfile.write(s.lower())
        s='#Columns '
        for key in mydata.data.keys():
            s=s+key+' '
        s=s+'\n'
        myoutfile.write(s)
        s=''
        for i in range(len(mydata.data[key])):
            for ckey in mydata.data:
                s=s+str(mydata.data[ckey][i])+' '
            s=s+'\n'
            myoutfile.write(s)
        myoutfile.close()
