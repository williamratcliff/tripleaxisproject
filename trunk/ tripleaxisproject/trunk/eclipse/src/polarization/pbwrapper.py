import ctypes
import cstruct
import numpy as N
# Defined ctypes
from ctypes import c_void_p, c_int, c_long, c_char, c_char_p,c_double
from ctypes import byref as _ref
c_void_pp = ctypes.POINTER(c_void_p)
c_char_pp = ctypes.POINTER(c_char_p)
c_int_p = ctypes.POINTER(c_int)
c_double_p=ctypes.POINTER(c_double)
c_ulong_p=ctypes.POINTER(c_ulong)


##class PBdatapt(ctypes.Structure):
##    _fields_=[("sec",c_long*4),
##            ("Y",c_double*4),
##            ("Yesq",c_double*4),
##            ("lambI",c_double),
##            ("lambF",c_double),
##            ("C",c_void_pp),
##            ("Cesq",c_void_pp),
##            ("S",c_double*4),
##            ("Sesq",c_double*4),
##            ("Nactive",c_int),
##            ("activeEq",c_int*4),
##            ("Nfree",c_int),
##            ("freeS",c_int*4)
##            ]


##PBdatapt = cstruct.cstruct(('sec',N.float64,(4,1)),('Y',N.float64,(4,1)),
##                    ('Yesq',N.float64,(4,1)),('lambI',N.float64,(4,1)),
##                    ('lamF',N.float64,(4,1)),('C',N.float64,(4,4)),
##                    ('Cesq',N.float64,(4,4)),('S',N.float64,(4,1)),
##                    ('Sesq',N.float64,(4,1)),('NActive',N.float64,(4,1)),
##                    ('activeEq',N.float64,(4,1)),('Nfree',N.intc),
##                    ('freeS',numpy.intc))


class PBindata(ctypes.Structure):
    _fields_=[("Ei",c_double_p),
            ("Ef",c_double_p),
            ("Cpp",c_double_p),
            ("Cmm",c_double_p),
            ("Cpm",c_double_p),
            ("Cmp",c_double_p),
            ("Epp",c_double_p),
            ("Emm",c_double_p),
            ("Epm",c_double_p),
            ("Emp",c_double_p),
            ("tpp",c_ulong_p),
            ("tmm",c_ulong_p),
            ("tpm",c_ulong_p),
            ("tmp",c_ulong_p)
            ]

#typedef struct {
#  double *Ei, *Ef ;
#  double *Cpp, *Cmm, *Cpm, *Cmp ;
#  double *Epp, *Emm, *Epm, *Emp ;
#  unsigned long *tpp, *tmm, *tpm, *tmp ;
#} PBindata ;

class PBoutdata(ctypes.Structure):
    _fields_=[("Spp",c_double_p),
            ("Smm",c_double_p),
            ("Spm",c_double_p),
            ("Smp",c_double_p),
            ("Epp",c_double_p),
            ("Emm",c_double_p),
            ("Epm",c_double_p),
            ("Emp",c_double_p)
            ]

#typedef struct {
#  double *Spp, *Smm, *Spm, *Smp ;
#  double *Epp, *Emm, *Epm, *Emp ;
#} PBoutdata ;

class PBflags(ctypes.Structure):
    _fields_=[("MonitorCorrect",c_int),
            ("PolMonitorCorrect",c_int),
            ("Debug",c_int),
            ("SimFlux",c_int),
            ("CountsEnable",c_int*4),
            ("CountsAdd1",c_int*4),
            ("CountsAdd2",c_int*4),
            ("Sconstrain",c_int*4),
            ("Spp",c_double*4),
            ("Smm",c_double*4),
            ("Spm",c_double*4),
            ("Smp",c_double*4)
            ]


#typedef struct {
#  int MonitorCorrect ;
#  int PolMonitorCorrect ;
#  int Debug ;
#  int SimFlux ;
#  int SimDeviate ;
#  int CountsEnable[4] ;
#  int CountsAdd1[4] ;
#  int CountsAdd2[4] ;
#  int Sconstrain[4] ;
#  double Spp[4], Smm[4], Spm[4], Smp[4] ;
#} PBflags ;


#/* entrypoints */

#int PBcorrectData(char *PCellFile, char *ACellFile, PBflags *flgs,
#		  int npts, PBindata *in, PBoutdata *out) ;
#int PBsim(char *filename) ;
#int PBreadflags(char *filename) ;


cell1="PcellBT7Jan072008.txt"
cell2="AcellBT7Jan72008.txt"

lib=ctypes.cdll['polarization.dll']
polcor=lib.PBcorrectData
flag=polcor(,ctypes.byref(argv))
