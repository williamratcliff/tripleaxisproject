import scriptutil as SU
import d10reader2 as d10reader

class DataSet(object):
    def __init__(self):
        self.data=[]
        self.hmin=[]
        self.hmax=[]   
        self.kmin=[]
        self.kmax=[]
        self.lmin=[]
        self.lmax=[]
        self.tth=[]
        self.omega_min=[]
        self.omega_max=[]
        self.tth_min=[]
        self.tth_max=[]
        self.filenums=[]

def read_files(mydirectory,myfilebase,myend):
    myfilebaseglob=myfilebase+'*.'+myend
    print myfilebaseglob
    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    #SU.printr(flist)
    dataset=DataSet()
    i=0
    for myfile in flist:
        print myfile
        data=d10reader.reader2(myfile)
        dataset.data.append(data)
        dataset.hmin.append(data.Hmin)
        dataset.kmin.append(data.Kmin)
        dataset.lmin.append(data.Lmin)
        dataset.hmax.append(data.Hmax)
        dataset.kmax.append(data.Kmax)
        dataset.lmax.append(data.Lmax)
        dataset.tth.append(data.tth)
        dataset.filenums.append(data.file_number)
        dataset.omega_min.append(min(data.angle1))
        dataset.omega_max.append(max(data.angle1))
        if i==0:
            om_min=min(data.angle1)
            om_max=max(data.angle1)
        else:
            om_min=min(om_min,min(data.angle1))
            om_max=max(om_max,max(data.angle1))
        
        try:
            dataset.tth_min.append(min(data.angle2))
        except:
            print 'oops min'
        try:
            dataset.tth_max.append(max(data.angle2))
        except:
            print 'oops max'
        i=i+1
    dataset.omega_minimum=om_min
    dataset.omega_maximum=om_max   
    return dataset
    
    
if __name__=="__main__":
    myfilebase=str(19262)
    myend='dat'
    mydirectory=r'c:\tbmno3\aug25_2009_ill'
    dataset=read_files(mydirectory,myfilebase,myend)
    print min(dataset.hmin)
    print 'done'
    
    