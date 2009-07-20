import numpy as N
import uncertainty
eps=1e-8

"""
Current notes:
My current thoughts are to have one giant table so that we can do sorting
easily.



"""

class Component(object):
        """This is the Component class.  A Component must have a name, for example, 'a1'
        Furthermore, it is given a set of values and stderr for initialization.
        units are optional.  Internally, we store a "measurement".  This can be
        accesed from measurement.x, measurement.dx
        """

        def _getx(self): return self.measurement.x
        def _setx(self,x): 
                self.measurement.x=x
        def _get_variance(self): return self.measurement.variance
        def _set_variance(self,variance):
                self.measurement.variance=variance
        def _getdx(self): return numpy.sqrt(self.variance)
        def _setdx(self,dx):
                # Direct operation
                #    variance = dx**2
                # Indirect operation to avoid temporaries
                self.variance[:] = dx
                self.variance **= 2
        x=property(_getx,_setx,doc='value')
        variance=property(_get_variance,_set_variance,doc='variance')
        dx = property(_getdx,_setdx,doc="standard deviation")

        #Out of laziness, I am defining properties of x, variance, and dx, but these only effect
        #the measurement objects attributes.  However, people should do operations on the motor object
        #NOT on these components, otherwise errors will not propagate correctly!!!
        def __init__(self, name,values,err, units=None,aliases=None):
                self.name=name
                self.aliases=aliases
                self.units=units
                #self.values=data
                #self.err=err
                self.measurement=uncertainty.Measurement(values,err**2)
        # Numpy array slicing operations
        def __len__(self):
                return len(self.x)
        def __getitem__(self,key):
                return Measurement(self.x[key],self.variance[key])
        def __setitem__(self,key,value):
                self.x[key] = value.x
                self.variance[key] = value.variance
        def __delitem__(self, key):
                del self.x[key]
                del self.variance[key]
        #def __iter__(self): pass # Not sure we need iter

        # Normal operations: may be of mixed type
        def __add__(self, other):
                if isinstance(other,Measurement):
                        return Measurement(*err1d.add(self.x,self.variance,other.x,other.variance))
                else:
                        return Measurement(self.x+other, self.variance+0) # Force copy
        def __sub__(self, other):
                if isinstance(other,Measurement):
                        return Measurement(*err1d.sub(self.x,self.variance,other.x,other.variance))
                else:
                        return Measurement(self.x-other, self.variance+0) # Force copy
        def __mul__(self, other):
                if isinstance(other,Measurement):
                        return Measurement(*err1d.mul(self.x,self.variance,other.x,other.variance))
                else:
                        return Measurement(self.x*other, self.variance*other**2)
        def __truediv__(self, other):
                if isinstance(other,Measurement):
                        return Measurement(*err1d.div(self.x,self.variance,other.x,other.variance))
                else:
                        return Measurement(self.x/other, self.variance/other**2)
        def __pow__(self, other):
                if isinstance(other,Measurement):
                        # Haven't calcuated variance in (a+/-da) ** (b+/-db)
                        return NotImplemented
                else:
                        return Measurement(*err1d.pow(self.x,self.variance,other))

        # Reverse operations
        def __radd__(self, other):
                return Measurement(self.x+other, self.variance+0) # Force copy
        def __rsub__(self, other):
                return Measurement(other-self.x, self.variance+0)
        def __rmul__(self, other):
                return Measurement(self.x*other, self.variance*other**2)
        def __rtruediv__(self, other):
                x,variance = err1d.pow(self.x,self.variance,-1)
                return Measurement(x*other,variance*other**2)
        def __rpow__(self, other): return NotImplemented

        # In-place operations: may be of mixed type
        def __iadd__(self, other):
                if isinstance(other,Measurement):
                        self.x,self.variance \
                            = err1d.add_inplace(self.x,self.variance,other.x,other.variance)
                else:
                        self.x+=other
                return self
        def __isub__(self, other):
                if isinstance(other,Measurement):
                        self.x,self.variance \
                            = err1d.sub_inplace(self.x,self.variance,other.x,other.variance)
                else:
                        self.x-=other
                return self
        def __imul__(self, other):
                if isinstance(other,Measurement):
                        self.x, self.variance \
                            = err1d.mul_inplace(self.x,self.variance,other.x,other.variance)
                else:
                        self.x *= other
                        self.variance *= other**2
                return self
        def __itruediv__(self, other):
                if isinstance(other,Measurement):
                        self.x,self.variance \
                            = err1d.div_inplace(self.x,self.variance,other.x,other.variance)
                else:
                        self.x /= other
                        self.variance /= other**2
                return self
        def __ipow__(self, other):
                if isinstance(other,Measurement):
                        # Haven't calcuated variance in (a+/-da) ** (b+/-db)
                        return NotImplemented
                else:
                        self.x,self.variance = err1d.pow_inplace(self.x, self.variance, other)
                return self

        # Use true division instead of integer division
        def __div__(self, other): return self.__truediv__(other)
        def __rdiv__(self, other): return self.__rtruediv__(other)
        def __idiv__(self, other): return self.__itruediv__(other)


        # Unary ops
        def __neg__(self):
                return Measurement(-self.x,self.variance)
        def __pos__(self):
                return self
        def __abs__(self):
                return Measurement(numpy.abs(self.x),self.variance)

        def __str__(self):
                #return str(self.x)+" +/- "+str(numpy.sqrt(self.variance))
                if numpy.isscalar(self.x):
                        return format_uncertainty(self.x,numpy.sqrt(self.variance))
                else:
                        return [format_uncertainty(v,dv)
                                for v,dv in zip(self.x,numpy.sqrt(self.variance))]
        def __repr__(self):
                return "Measurement(%s,%s)"%(str(self.x),str(self.variance))

        # Not implemented
        def __floordiv__(self, other): return NotImplemented
        def __mod__(self, other): return NotImplemented
        def __divmod__(self, other): return NotImplemented
        def __mod__(self, other): return NotImplemented
        def __lshift__(self, other): return NotImplemented
        def __rshift__(self, other): return NotImplemented
        def __and__(self, other): return NotImplemented
        def __xor__(self, other): return NotImplemented
        def __or__(self, other): return NotImplemented

        def __rfloordiv__(self, other): return NotImplemented
        def __rmod__(self, other): return NotImplemented
        def __rdivmod__(self, other): return NotImplemented
        def __rmod__(self, other): return NotImplemented
        def __rlshift__(self, other): return NotImplemented
        def __rrshift__(self, other): return NotImplemented
        def __rand__(self, other): return NotImplemented
        def __rxor__(self, other): return NotImplemented
        def __ror__(self, other): return NotImplemented

        def __ifloordiv__(self, other): return NotImplemented
        def __imod__(self, other): return NotImplemented
        def __idivmod__(self, other): return NotImplemented
        def __imod__(self, other): return NotImplemented
        def __ilshift__(self, other): return NotImplemented
        def __irshift__(self, other): return NotImplemented
        def __iand__(self, other): return NotImplemented
        def __ixor__(self, other): return NotImplemented
        def __ior__(self, other): return NotImplemented

        def __invert__(self): return NotImplmented  # For ~x
        def __complex__(self): return NotImplmented
        def __int__(self): return NotImplmented
        def __long__(self): return NotImplmented
        def __float__(self): return NotImplmented
        def __oct__(self): return NotImplmented
        def __hex__(self): return NotImplmented
        def __index__(self): return NotImplmented
        def __coerce__(self): return NotImplmented

        def log(self):
                return Measurement(*err1d.log(self.x,self.variance))

        def exp(self):
                return Measurement(*err1d.exp(self.x,self.variance))

        def log(val): return self.log()
        def exp(val): return self.exp()


class Motor(Component):
        """This is the motor class.  A Motor must have a name, for example, 'a1'
        Furthermore, it is given a set of values and stderr for initialization.
        units are optional.  Internally, we store a "measurement".  This can be
        accesed from measurement.x, measurement.dx
        """


        def __init__(self,name,values=None,err=None,units='degrees',isDistinct=True, 
                     window=eps,aliases=None,friends=None,spectator=False
                     isInterpolatable=False):
                self.name=name
                self.units=units
                self.measurement=uncertainty.Measurement(values, err**2)
                self.aliases=aliases
                self.isDistinct=isDistinct
                self.isInterpolatable=isInterpolatable
                self.spectator=spectator
                #The spectator flag says if I was moving or not during a scan.
                self.window=window #The window in which values are to be deemed equal, it assumes that the isDistinct flag is set 
                self.friends=friends  
                #If I am updated, then my friends might need to be updated, for example, hkl->


class SampleEnvironment(Component):
        """This is the motor class.  A Motor must have a name, for example, 'a1'
        Furthermore, it is given a set of values and stderr for initialization.
        units are optional.  Internally, we store a "measurement".  This can be
        accesed from measurement.x, measurement.dx
        """


        def __init__(self,name,values=None,err=None,units='degrees',isDistinct=True, 
                     window=eps,aliases=None,friends=None,spectator=False
                     isInterpolatable=False):
                self.name=name
                self.units=units
                self.measurement=uncertainty.Measurement(values, err**2)
                self.aliases=aliases
                self.isDistinct=isDistinct
                self.isInterpolatable=isInterpolatable
                self.spectator=spectator
                #The spectator flag says if I was moving or not during a scan.
                self.window=window #The window in which values are to be deemed equal, it assumes that the isDistinct flag is set 
                self.friends=friends  
                #If I am updated, then my friends might need to be updated, for example, hkl->




class Mosaic(object):
        def __init__(self,horizontal, vertical=None):
                self.horizontal=horizontal
                self.vertical=vertical

class Monochromator(object):
        """A monochromator"""
        def __init__(self, name='Monochromator', 
                     vertical_focus=None,
                     horizontal_focus=None,
                     blades=None,
                     mosaic=None,
                     aliases=None,
                     dspacing=None  #dspacing of the monochromator
                     ):
                self.name=name
                self.vertical_focus=vertical_focus
                self.horizontal_focus=horizontal_focus
                sellf.blades=blades
                self.mosaic=mosaic
                self.focus_cu=Motor('focus_cu',values=None,err=None,units='degrees',isDistinct=True,
                                    isInterpolatable=True)
                self.focus_pg=Motor('focus_pg',values=None,err=None,units='degrees',isDistinct=True,
                                    isInterpolatable=True)

class Filters(object):
        """Filters"""
        def __init__(self):
                self.filter_rotation=Motor('filter_rotation',values=None,err=None,units='degrees',isDistinct=False,
                                           isInterpolatable=False)
                self.filter_tilt=Motor('filter_tilt',values=None,err=None,units='degrees',isDistinct=False,
                                       isInterpolatable=False)
                self.filter_translation=Motor('filter_translation',values=None,err=None,units='none',isDistinct=True,
                                              isInterpolatable=False)
                #filter translation has values of "IN" and "OUT" 

class Apertures(object):
        def __init__(self):
                self.aperture_horizontal=Motor('aperture_horizontal',values=None,err=None,units='degrees',isDistinct=False,
                                               isInterpolatable=False)
                self.aperture_vertical=Motor('aperture_vertical',values=None,err=None,units='degrees',isDistinct=False,
                                             isInterpolatable=False)

class Time(object):
        def __init__(self):
                self.duration=Motor('duration',values=None,err=None,units='seconds',isDistinct=True,
                                                       isInterpolatable=True)
                self.timestamp=Motor('timestamp',values=None,err=None,units='seconds',isDistinct=False,
                                                       isInterpolatable=True)
                #We should decide if we ignore timestamps, for now let's do so....

class Temperature(object):                
        def __init__(self):
                self.temperature=SampleEnvironment('temperature',values=None,err=None,units='K',isDistinct=True,
                                               isInterpolatable=True)
                self.temperature_control_reading=SampleEnvironment('temperature_control_reading',values=None,err=None,units='K',isDistinct=False,
                                               isInterpolatable=True)
                self.temperature_heater_power=SampleEnvironment('temperature_heater_power',values=None,err=None,units='percentage',isDistinct=False,
                                               isInterpolatable=True) #usuallly a percentage for most controllers...
                self.temperature_setpoint=SampleEnvironment('temperature_setpoint',values=None,err=None,units='K',isDistinct=False,
                                               isInterpolatable=True)
                

        
                
                
class Analyzer(object):
        """An analyzer"""
        def __init__(self, name='Analyzer',
                     vertical_focus=None,
                     horizontal_focus=None,
                     blades=None,
                     mosaic=None,
                     dspacing=None  #dspacing of the monochromator
                     ):
                self.name=name
                self.blades=blades
                self.mosaic=mosaic




class Primary_Motors(object):
        def __init__(self):
                self.a1=Motor('a1',values=None,err=None,units='degrees',isDistinct=True, isInterpolatable=True)
                self.a2=Motor('a2',values=None,err=None,units='degrees',isDistinct=True, isInterpolatable=True)
                self.a3=Motor('a3',values=None,err=None,units='degrees',isDistinct=True, isInterpolatable=True)
                self.a4=Motor('a4',values=None,err=None,units='degrees',isDistinct=True, isInterpolatable=True)
                self.a5=Motor('a5',values=None,err=None,units='degrees',isDistinct=True, isInterpolatable=True)
                self.a6=Motor('a6',values=None,err=None,units='degrees',isDistinct=True, isInterpolatable=True)
                self.sample_elevator=Motor('sample_elevation',values=None,err=None,units='degrees',isDistinct=False
                                           , isInterpolatable=True)
                self.sample_elevator=Motor('a6',values=None,err=None,units='degrees',isDistinct=False
                                           , isInterpolatable=True)
                self.dfm_rotation=Motor('dfm_rotation',values=None,err=None,units='degrees',isDistinct=True
                                        , isInterpolatable=True)
                self.analyzer_rotation=Motor('analyzer_rotation',values=None,err=None,units='degrees',isDistinct=True
                                             , isInterpolatable=True)
                self.aperture_horizontal=Motor('aperture_horizontal',values=None,err=None,units='degrees',isDistinct=True
                                               , isInterpolatable=True)
                self.aperture_vertical=Motor('aperture_vertical',values=None,err=None,units='degrees',isDistinct=True
                                             , isInterpolatable=True)


class Physical_Motors(object):
        def __init__(self):
                self.h=Motor('h',values=None,err=None,units='rlu',isDistinct=True,
                             isInterpolatable=True)
                self.k=Motor('k',values=None,err=None,units='rlu',isDistinct=True,
                             isInterpolatable=True)
                self.l=Motor('l',values=None,err=None,units='rlu',isDistinct=True,
                             isInterpolatable=True)
                self.e=Motor('e',values=None,err=None,units='rlu',isDistinct=True,
                             isInterpolatable=True)
                self.qx=Motor('qx',values=None,err=None,units='rlu',isDistinct=True,
                              isInterpolatable=True)
                self.qy=Motor('qy',values=None,err=None,units='rlu',isDistinct=True,
                              isInterpolatable=True)
                self.qz=Motor('qz',values=None,err=None,units='rlu',isDistinct=True,
                              isInterpolatable=True)
                self.hkl=Motor('hkl',values=None,err=None,units='rlu',isDistinct=True,
                               isInterpolatable=True)

class Collimators(object):
        """Our collimators:  
        post_analyzer_collimator:user_defined
        post_monochromator_collimator:user_defined
        pre_analyzer_collimator:user_defined
        pre_monochromator_collimator:measured
        For the radial and soller collimators, these are angles on the track.  Others
        self.radial_collimator:  For horizontal focusing mode, rc=a6
        self.soller_collimator:  When in flat mode, the 50 minute collimator=a6, 
        others have a fixed offset relative to this one.  Within some window, these
        are likely to remain constant.  Thus, from the angle on the track, we could probably
        determine what the likely collimation status is.  This may be a candidate for a 
        "slowly varying" property.  For now, we will treat these as nonDistict

        """
        def __init__(self):
                self.post_analyzer_collimator=Motor('post_analyzer_collimator',values=None,err=None,units='minutes',isDistinct=False)
                self.post_monochromator_collimator=Motor('post_monochromator_collimator',values=None,err=None,units='minutes',isDistinct=False)
                self.pre_analyzer_collimator=Motor('pre_analyzer_collimator',values=None,err=None,units='minutes',isDistinct=False)
                self.pre_monochromator_collimator=Motor('post_monochromator_collimator',values=None,err=None,units='minutes',isDistinct=True)
                self.radial_collimator=Motor('radial_collimator',values=None,err=None,units='degrees',isDistinct=True, window=2.0)
                self.soller_collimator=Motor('soller_collimator',values=None,err=None,units='degrees',isDistinct=False, window=2.0)



class TripleAxis(object):
        def __init__(self):
                self.monochromator=Monochromator()
                self.analyzer=Analyzer()



a1=Motor('a1',values=None,err=None,units='degrees',isdistinct=True)
a2=Motor('a2',values=None,err=None,units='degrees',isdistinct=True)
a3=Motor('a3',values=None,err=None,units='degrees',isdistinct=True)
a4=Motor('a4',values=None,err=None,units='degrees',isdistinct=True)
a5=Motor('a5',values=None,err=None,units='degrees',isdistinct=True)
a6=Motor('a6',values=None,err=None,units='degrees',isdistinct=True)
sample_elevator=Motor('sample_elevation',values=None,err=None,units='degrees',isdistinct=False)
sample_elevator=Motor('a6',values=None,err=None,units='degrees',isdistinct=False)
dfm_rotation=Motor('dfm_rotation',values=None,err=None,units='degrees',isdistinct=True)
analyzer_rotation=Motor('analyzer_rotation',values=None,err=None,units='degrees',isdistinct=True)
aperture_horizontal=Motor('aperture_horizontal',values=None,err=None,units='degrees',isdistinct=True)
aperture_vertical=Motor('aperture_vertical',values=None,err=None,units='degrees',isdistinct=True)

def data_abstraction_layer(self):
        self.metadata={}
        self.additional_metadata={}
        self.metadata['count_info']={}
        self.metadata['count_info']['monitor_base']=None #float(tokenized[6])
        self.metadata['count_info']['monitor_prefactor']=None#float(tokenized[7])
        self.metadata['count_info']['monitor']=None#self.metadata['count_info']['monitor_base']*self.metadata['count_info']['monitor_prefactor']
        self.metadata['count_info']['count_type']=None  #can be 'monitor', 'time' #tokenized[8].strip("'").lower()
        self.metadata['count_info']['signal']=None  #for example, 'detector'
        self.metadata['count_info']['varying']=None
        self.metadata['count_info']['ranges']=None
        self.metadata['count_info']['analyzerdetectormode']=None
        self.metadata['count_info']['AnalyzerDetectorDevicesOfInterest'.lower()]=None
        self.metadata['count_info']['AnalyzerDDGroup'.lower()]=None
        self.metadata['count_info']['AnalyzerPSDGroup'.lower()]=None
        self.metadata['count_info']['AnalyzerSDGroup'.lower()]=None
        self.metadata['count_info']['AnalyzerDoorDetectorGroup'.lower()]=None
        self.metadata['count_info']['analyzerfocusmode']=None
        self.metadata['count_info']['monovertifocus']=None
        self.metadata['count_info']['monohorizfocus']=None


        self.metadata['file_info']={}
        self.metadata['file_info']['filename']=None#tokenized[0].strip("'")
        self.metadata['file_info']['filebase']=None#self.metadata['file_info']['filename'][0:5]
        self.metadata['file_info']['fileseq_number']=None
        self.metadata['file_info']['scantype']=None#tokenized[5].strip("'").lower()
        self.metadata['file_info']['instrument']=None#self.metadata['file_info']['filename'].split('.')[1].lower()
        self.metadata['file_info']['comment']=None #myfile.readline().rstrip()
        self.metadata['file_info']['scan_description']=None
        self.metadata['file_info']['experiment_id']=None
        self.metadata['file_info']['fixed_devices']=None

        self.metadata['timestamp']={}
        self.metadata['timestamp']['month']=None#int, for icp data it is translated using the months dict
        self.metadata['timestamp']['day']=None#int
        self.metadata['timestamp']['year']=None#int
        self.metadata['timestamp']['time']=None#str
        self.metadata['timestamp']['epoch']=None#float

        self.metadata['collimations']={}
        self.metadata['collimations']['coll1']=None#float(tokenized[0])
        self.metadata['collimations']['coll2']=None#float(tokenized[1])
        self.metadata['collimations']['coll3']=None#float(tokenized[2])
        self.metadata['collimations']['coll4']=None#float(tokenized[3])

        self.metadata['mosaic']={}
        self.metadata['mosaic']['mosaic_monochromator']=None#float(tokenized[4])
        self.metadata['mosaic']['mosaic_sample']=None#float(tokenized[5])
        self.metadata['mosaic']['mosaic_analyzer']=None#float(tokenized[6])

        self.metadata['dspacing']={}
        self.metadata['dspacing']['monochromator_dspacing']=None#float(tokenized[3])
        self.metadata['dspacing']['analyzer_dspacing']=None#float(tokenized[4])


        self.metadata['energy_info']={}
        self.metadata['energy_info']['wavelength']=None#float(tokenized[7])
        self.metadata['energy_info']['ef']=None#float(tokenized[2])
        self.metadata['energy_info']['efixed']=None#Should be 'ei,ef'#tokenized[4]


        self.metadata['orient1']={}
        self.metadata['orient1']['h']=None#float(tokenized[7])
        self.metadata['orient1']['k']=None#float(tokenized[8])
        self.metadata['orient1']['l']=None#float(tokenized[9])
        #ignore "angle" field
        self.metadata['orient2']={}
        self.metadata['orient2']['h']=None#float(tokenized[11])
        self.metadata['orient2']['k']=None#float(tokenized[12])
        self.metadata['orient2']['l']=None#float(tokenized[13])

        self.metadata['lattice']={}
        self.metadata['lattice']['a']=None#float(tokenized[0])
        self.metadata['lattice']['b']=None#float(tokenized[1])
        self.metadata['lattice']['c']=None#float(tokenized[2])
        self.metadata['lattice']['alpha']=None#float(tokenized[3])
        self.metadata['lattice']['beta']=None#float(tokenized[4])
        self.metadata['lattice']['gamma']=None#float(tokenized[5])

        self.metadata['q_center']={}
        self.metadata['q_center']['e_center']=None#float(tokenized[0])
        self.metadata['q_center']['h_center']=None#(tokenized[0])
        self.metadata['q_center']['k_center']=None#(tokenized[1])
        self.metadata['q_center']['l_center']=None#(tokenized[2])


        self.metadata['q_step']={}
        self.metadata['q_step']['delta_h']=None#float(tokenized[3])
        self.metadata['q_step']['delta_k']=None#float(tokenized[4])
        self.metadata['q_step']['delta_l']=None#float(tokenized[5])
        self.metadata['q_step']['delta_e']=None#float(tokenized[1])


        self.metadata['temperature_info']={}
        self.metadata['temperature_info']['Tstart']=None#float(tokenized[5])
        self.metadata['temperature_info']['Tstep']=None#float(tokenized[6])
        self.metadata['temperature_info']['units']=None#float(tokenized[6])

        self.metadata['magnetic_field']={}
        self.metadata['magnetic_field']['hfield']=None#float(tokenized[10])
        return






