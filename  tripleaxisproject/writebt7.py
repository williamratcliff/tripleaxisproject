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
        myoutfile=open(myoutfilestr,'wt')
        for i in range(len(mydata.header)):
            s=mydata.header[i]+'\n'
            myoutfile.write(s)
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
