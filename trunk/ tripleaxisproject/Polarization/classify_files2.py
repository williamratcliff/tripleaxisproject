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
        #self.pol_info=pol_info()
        #self.pm=pol_info()
        #self.mp=pol_info()
        #self.mm=pol_info()
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
        data['full_data']=mydata
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
                data['temp']={}
                data['temp']['min']=temp.min()
                data['temp']['max']=temp.max()
                data['temp']['center']=N.average(temp)
            except ValueError:
                pass
        if mydata.data.has_key('magfield'):
            try:
                magfield=N.array(mydata.data['magfield'],'float64')
                #data['environment']['has_temp']=True
                data['magfield']={}
                data['magfield']['min']=magfield.min()
                data['magfield']['max']=magfield.max()
                data['magfield']['center']=N.average(magfield)
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
        data['polarization state']=key
        mycatalog.data.append(data)
        mycatalog.files.append(data['filename'])
        curr_catalog=mycatalog
        if data.has_key('h'):
            if curr_catalog.h_range.max==None:
                curr_catalog.h_range.max=data['h']['max']
            else:
                curr_catalog.h_range.max=max(data['h']['max'],curr_catalog.h_range.max)
            if curr_catalog.h_range.min==None:
                curr_catalog.h_range.min=data['h']['min']
            else:
                curr_catalog.h_range.min=min(data['h']['min'],curr_catalog.h_range.min)
        if data.has_key('k'):
            if curr_catalog.k_range.max==None:
                curr_catalog.k_range.max=data['k']['max']
            else:
                curr_catalog.k_range.max=max(data['k']['max'],curr_catalog.k_range.max)
            if curr_catalog.k_range.min==None:
                curr_catalog.k_range.min=data['k']['min']
            else:
                curr_catalog.k_range.min=min(data['k']['min'],curr_catalog.k_range.min)
        if data.has_key('l'):
            if curr_catalog.l_range.max==None:
                curr_catalog.l_range.max=data['l']['max']
            else:
                curr_catalog.l_range.max=max(data['l']['max'],curr_catalog.l_range.max)
            if curr_catalog.l_range.min==None:
                curr_catalog.l_range.min=data['l']['min']
            else:
                curr_catalog.l_range.min=min(data['l']['min'],curr_catalog.l_range.min)

        if data.has_key('e'):
            if curr_catalog.e_range.max==None:
                curr_catalog.e_range.max=data['e']['max']
            else:
                curr_catalog.e_range.max=max(data['e']['max'],curr_catalog.e_range.max)
            if curr_catalog.e_range.min==None:
                curr_catalog.e_range.min=data['e']['min']
            else:
                curr_catalog.e_range.min=min(data['e']['min'],curr_catalog.e_range.min)

        if data.has_key('a3'):
            if curr_catalog.a3_range.max==None:
                curr_catalog.a3_range.max=data['a3']['max']
            else:
                curr_catalog.a3_range.max=max(data['a3']['max'],curr_catalog.a3_range.max)
            if curr_catalog.a3_range.min==None:
                curr_catalog.a3_range.min=data['a3']['min']
            else:
                curr_catalog.a3_range.min=min(data['a3']['min'],curr_catalog.a3_range.min)


        if data.has_key('a4'):
            if curr_catalog.a4_range.max==None:
                curr_catalog.a4_range.max=data['a4']['max']
            else:
                curr_catalog.a4_range.max=max(data['a4']['max'],curr_catalog.a4_range.max)
            if curr_catalog.a4_range.min==None:
                curr_catalog.a4_range.min=data['a4']['min']
            else:
                curr_catalog.a4_range.min=min(data['a4']['min'],curr_catalog.a4_range.min)


        if data.has_key('temp'):
            if curr_catalog.temp_range.max==None:
                curr_catalog.temp_range.max=data['temp']['max']
            else:
                curr_catalog.temp_range.max=max(data['temp']['max'],curr_catalog.temp_range.max)
            if curr_catalog.temp_range.min==None:
                curr_catalog.temp_range.min=data['temp']['min']
            else:
                curr_catalog.temp_range.min=min(data['temp']['min'],curr_catalog.temp_range.min)

        if data.has_key('magfield'):
            if curr_catalog.magfield_range.max==None:
                curr_catalog.magfield_range.max=data['magfield']['max']
            else:
                curr_catalog.magfield_range.max=max(data['magfield']['max'],curr_catalog.magfield_range.max)
            if curr_catalog.magfield_range.min==None:
                curr_catalog.magfield_range.min=data['magfield']['min']
            else:
                curr_catalog.magfield_range.min=min(data['magfield']['min'],curr_catalog.magfield_range.min)




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
    #print mycatalog.pm.files[0]
    #print mycatalog.pm.data[0]['h']
    #print mycatalog.pm.data[0]['k']
    #print mycatalog.pm.data[0]['l']
    #print mycatalog.pm.data[0]['e']
    #print mycatalog.pm.data[0]['ei']
    #print mycatalog.pm.data[0]['a3']
    #print mycatalog.pm.data[0]['a4']
    #print mycatalog.pm.data[0]['environment']
    #print mycatalog.pm.data[0]['count_type']
#    print mycatalog.pm.temp_range.min
    print mycatalog.h_range.max
    print mycatalog.h_range.min




