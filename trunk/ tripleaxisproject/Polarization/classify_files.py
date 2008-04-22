import numpy as N
import scriptutil as SU
import re
import readncnr2 as readncnr
import sys,os
import math
threshold=1.0e-1


def autovectorized(f):
     """Function decorator to do vectorization only as necessary.
     vectorized functions fail for scalar inputs."""
     def wrapper(input):
         if N.isscalar(input)==False:
             return N.vectorize(f)(input)
         return f(input)
     return wrapper



@autovectorized
def myradians(x):
    return math.radians(x)

def calc_energy(angle,dspace):
    anglerad=myradians(angle)
    tau=2*N.pi/dspace
    k=tau/N.sqrt(2-2*N.cos(anglerad))
    energy=2.072142*k*k
    return energy



#plus refers to a flipper state that is on
#minus refers to a flipper state that is off


class field_range():
    def __init__(self):
        self.min=None
        self.max=None



class pol_info:
    def __init__(self):
        self.data=[]
        self.files=[]
        self.h_range=field_range()
        self.k_range=field_range()
        self.l_range=field_range()
        self.e_range=field_range()
        self.a3_range=field_range()
        self.a4_range=field_range()
        self.temp_range=field_range()
        self.magfield_range=field_range()
        return


class catalog:
    def __init__(self):
        self.pp=pol_info()
        self.pm=pol_info()
        self.mp=pol_info()
        self.mm=pol_info()
        return

def readfiles(flist):
    #SU.printr(flist)
    mydatareader=readncnr.datareader()
    mycatalog=catalog()
    half_polarized=0
    for currfile in flist:
        #print currfile
        key_i='m'
        hsample='none'
        vsample='none'
        mydata=mydatareader.readbuffer(currfile)
        #mydata=mydatareader.readbuffer(currfile,lines=3)
        if mydata.data.has_key('hsample'):
            hsample=mydata.data['hsample'][-1]
        if mydata.data.has_key('vsample'):
            vsample=mydata.data['vsample'][-1]
        ei_exists=0
        if mydata.data.has_key('eiflip'):
            ei_exists=1
            #print mydata.data['eiflip']
            if N.average(mydata.data['eiflip'])>threshold:
                key_i='p'
        key_f='m'
        ef_exists=0
        if mydata.data.has_key('efflip'):
            ef_exists=1
            if N.average(mydata.data['efflip'])>threshold:
                key_f='p'
        key=key_i+key_f
        fully_polarized=ei_exists*ef_exists  #is 1 if both cells are present
        data={}
        data['absolute_filename']=currfile
        data['fully_polarized']=fully_polarized
        data['hsample']=hsample
        data['vsample']=vsample

        if mydata.data.has_key('qx'):
            try:
                h=N.array(mydata.data['qx'],'float64')
                data['h']={}
                data['h']['min']=h.min()
                data['h']['max']=h.max()
                data['h']['center']=N.average(h)
            except ValueError:
                pass
                #data['has_h']=False
        if mydata.data.has_key('qy'):
            try:
                k=N.array(mydata.data['qy'],'float64')
                data['k']={}
                data['k']['min']=k.min()
                data['k']['max']=k.max()
                data['k']['center']=N.average(k)
            except ValueError:
                pass
        if mydata.data.has_key('qz'):
            try:
                l=N.array(mydata.data['qz'],'float64')
                data['l']={}
                data['l']['min']=l.min()
                data['l']['max']=l.max()
                data['l']['center']=N.average(l)
            except ValueError:
                pass
                #data['has_l']=False
        if mydata.data.has_key('e'):
            try:
                e=N.array(mydata.data['e'],'float64')
                data['e']={}
                data['e']['min']=e.min()
                data['e']['max']=e.max()
                data['e']['center']=N.average(e)
            except ValueError:
                pass
                #data['has_e']=False
        data['environment']={}
        if mydata.data.has_key('temp'):
            try:
                temp=N.array(mydata.data['temp'],'float64')
                #data['environment']['has_temp']=True
                data['environment']['temp']={}
                data['environment']['temp']['min']=temp.min()
                data['environment']['temp']['max']=temp.max()
                data['environment']['temp']['center']=N.average(temp)
            except ValueError:
                pass
        if mydata.data.has_key('magfield'):
            try:
                magfield=N.array(mydata.data['magfield'],'float64')
                #data['environment']['has_temp']=True
                data['environment']['magfield']={}
                data['environment']['magfield']['min']=magfield.min()
                data['environment']['magfield']['max']=magfield.max()
                data['environment']['magfield']['center']=N.average(magfield)
            except ValueError:
                pass

        if mydata.data.has_key('a2'):
            a2=N.array(mydata.data['a2'])
            dspace_m=mydata.metadata['dspacing']['monochromator_dspacing']
            ei=calc_energy(a2,dspace_m)
            data['ei']={}
            data['ei']['min']=ei.min()
            data['ei']['max']=ei.max()
            data['ei']['center']=N.average(ei)

        if mydata.data.has_key('a3'):
            try:
                a3=N.array(mydata.data['a3'],'float64')
                #data['environment']['has_temp']=True
                data['a3']={}
                data['a3']['min']=a3.min()
                data['a3']['max']=a3.max()
                data['a3']['center']=N.average(a3)
            except ValueError:
                pass

        if mydata.data.has_key('a4'):
            try:
                a4=N.array(mydata.data['a4'],'float64')
                #data['environment']['has_temp']=True
                data['a4']={}
                data['a4']['min']=a4.min()
                data['a4']['max']=a4.max()
                data['a4']['center']=N.average(a4)
            except ValueError:
                pass

        try:
            data['count_type']=mydata.metadata['count_info']['count_type']
        except KeyError:
            pass

        pattern = re.compile('^(?P<base>[^.]*?)(?P<seq>[0-9]*)(?P<ext>[.].*)?$')
        filename=os.path.split(currfile)[-1]
        #print filename
        match = pattern.match(filename)
        dict((a,match.group(a)+"") for a in ['base','seq','ext'])
        data['filebase']=match.group('base')
        data['fileseq_number']=match.group('seq')
        data['filename']=data['filebase']+data['fileseq_number']

        if key=='pp':
            mycatalog.pp.data.append(data)
            mycatalog.pp.files.append(data['filename'])
            if data.has_key('h'):
                mycatalog.pp.h_range.max=max(mycatalog.pp.h_range.max,data['h']['max'])
                mycatalog.pp.h_range.min=min(mycatalog.pp.h_range.min,data['h']['min'])
            if data.has_key('k'):
                mycatalog.pp.k_range.max=max(mycatalog.pp.k_range.max,data['k']['max'])
                mycatalog.pp.k_range.min=min(mycatalog.pp.k_range.min,data['k']['min'])
            if data.has_key('l'):
                mycatalog.pp.l_range.max=max(mycatalog.pp.l_range.max,data['l']['max'])
                mycatalog.pp.l_range.min=min(mycatalog.pp.l_range.min,data['l']['min'])
            if data.has_key('e'):
                mycatalog.pp.e_range.max=max(mycatalog.pp.e_range.max,data['e']['max'])
                mycatalog.pp.e_range.min=min(mycatalog.pp.e_range.min,data['e']['min'])
            if data.has_key('temp'):
                mycatalog.pp.temp_range.max=max(mycatalog.pp.temp_range.max,data['temp']['max'])
                mycatalog.pp.temp_range.min=min(mycatalog.pp.temp_range.min,data['temp']['min'])
            if data.has_key('magfield'):
                mycatalog.pp.magfield_range.max=max(mycatalog.pp.magfield_range.max,data['magfield']['max'])
                mycatalog.pp.magfield_range.min=min(mycatalog.pp.magfield_range.min,data['magfield']['min'])


        elif key=='pm':
            mycatalog.pm.data.append(data)
            mycatalog.pm.files.append(data['filename'])
            if data.has_key('h'):
                curr_catalog=mycatalog.pm
                if curr_catalog.h_range.max==None:
                    curr_catalog.h_range.max=data['h']['max']
                else:
                    curr_catalog.h_range.max=max(data['h']['max'],curr_catalog.h_range.max)
                if curr_catalog.h_range.min==None:
                    curr_catalog.h_range.min=data['h']['min']
                else:
                    curr_catalog.h_range.min=min(data['h']['min'],curr_catalog.h_range.min)


        elif key=='mp':
            mycatalog.mp.data.append(data)
            mycatalog.mp.files.append(data['filename'])
            if data.has_key('h'):
                mycatalog.mp.h_range.max=max(mycatalog.mp.h_range.max,data['h']['max'])
                mycatalog.mp.h_range.min=min(mycatalog.mp.h_range.min,data['h']['min'])
            if data.has_key('k'):
                mycatalog.mp.k_range.max=max(mycatalog.mp.k_range.max,data['k']['max'])
                mycatalog.mp.k_range.min=min(mycatalog.mp.k_range.min,data['k']['min'])
            if data.has_key('l'):
                mycatalog.mp.l_range.max=max(mycatalog.mp.l_range.max,data['l']['max'])
                mycatalog.mp.l_range.min=min(mycatalog.mp.l_range.min,data['l']['min'])
            if data.has_key('e'):
                mycatalog.mp.e_range.max=max(mycatalog.mp.e_range.max,data['e']['max'])
                mycatalog.mp.e_range.min=min(mycatalog.mp.e_range.min,data['e']['min'])
            if data.has_key('temp'):
                mycatalog.mp.temp_range.max=max(mycatalog.mp.temp_range.max,data['temp']['max'])
                mycatalog.mp.temp_range.min=min(mycatalog.mp.temp_range.min,data['temp']['min'])
            if data.has_key('magfield'):
                mycatalog.mp.magfield_range.max=max(mycatalog.mp.magfield_range.max,data['magfield']['max'])
                mycatalog.mp.magfield_range.min=min(mycatalog.mp.magfield_range.min,data['magfield']['min'])

        elif key=='mm':
            mycatalog.mm.data.append(data)
            mycatalog.mm.files.append(data['filename'])
            if data.has_key('h'):
                mycatalog.mm.h_range.max=max(mycatalog.mm.h_range.max,data['h']['max'])
                mycatalog.mm.h_range.min=min(mycatalog.mm.h_range.min,data['h']['min'])
            if data.has_key('k'):
                mycatalog.mm.k_range.max=max(mycatalog.mm.k_range.max,data['k']['max'])
                mycatalog.mm.k_range.min=min(mycatalog.mm.k_range.min,data['k']['min'])
            if data.has_key('l'):
                mycatalog.mm.l_range.max=max(mycatalog.mm.l_range.max,data['l']['max'])
                mycatalog.mm.l_range.min=min(mycatalog.mm.l_range.min,data['l']['min'])
            if data.has_key('e'):
                mycatalog.mm.e_range.max=max(mycatalog.mm.e_range.max,data['e']['max'])
                mycatalog.mm.e_range.min=min(mycatalog.mm.e_range.min,data['e']['min'])
            if data.has_key('temp'):
                mycatalog.mm.temp_range.max=max(mycatalog.mm.temp_range.max,data['temp']['max'])
                mycatalog.mm.temp_range.min=min(mycatalog.mm.temp_range.min,data['temp']['min'])
            if data.has_key('magfield'):
                mycatalog.mm.magfield_range.max=max(mycatalog.mm.magfield_range.max,data['magfield']['max'])
                mycatalog.mm.magfield_range.min=min(mycatalog.mm.magfield_range.min,data['magfield']['min'])

        #print currfile, key
    return mycatalog



if __name__=='__main__':
    myend='bt7'
    mydirectory=r'c:\bifeo3xtal\jan8_2008\9175'
    myfilebase=''
    myfilebaseglob=myfilebase+'*.'+myend
    #print myfilebaseglob
    flist = SU.ffind(mydirectory, shellglobs=(myfilebaseglob,))
    mycatalog=readfiles(flist)
    print mycatalog.pm.files[0]
    print mycatalog.pm.data[0]['h']
    print mycatalog.pm.data[0]['k']
    print mycatalog.pm.data[0]['l']
    print mycatalog.pm.data[0]['e']
    print mycatalog.pm.data[0]['ei']
    print mycatalog.pm.data[0]['a3']
    print mycatalog.pm.data[0]['a4']
    print mycatalog.pm.data[0]['environment']
    print mycatalog.pm.data[0]['count_type']
#    print mycatalog.pm.temp_range.min
    print mycatalog.pm.h_range.max
    print mycatalog.pm.h_range.min




