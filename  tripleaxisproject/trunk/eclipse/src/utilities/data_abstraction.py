import numpy as N


class Component(object):
        def __init__(self, name=None,units=None, data=None, err=None):
                self.name=name
                self.units=units
                self.data=data
                self.err=err
                
class Motor(Component):
        def __init__(self,name=None,units='degrees',data=None,err=None):
                super(Motor,self).__init__(name,units,data,err)
                

class Monitor(Component):
        """A monitor
        def __init__(self, name='Monitor1',
                     signaltype=None,
                     vertical_focus=None,
                     horizontal_focus=None,
                     Blades=None,
                     mosaic_vertical=None,
                     mosaic_horizntal=None,
                     dspacing=None  #dspacing of the monochromator
                     )
                




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
