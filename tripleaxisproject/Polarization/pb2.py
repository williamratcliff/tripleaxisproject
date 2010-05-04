#int PBcorrect(PBdatapt *d)
import ctypes
import cstruct
import numpy as N
# Defined ctypes
from ctypes import c_void_p, c_int, c_long, c_char, c_char_p,c_double
from ctypes import byref as _ref
c_void_pp = ctypes.POINTER(c_void_p)
c_char_pp = ctypes.POINTER(c_char_p)
c_int_p = ctypes.POINTER(c_int)


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


PBdatapt = cstruct.cstruct(('sec',N.float64,(4,1)),('Y',N.float64,(4,1)),
                    ('Yesq',N.float64,(4,1)),('lambI',N.float64,(4,1)),
                    ('lamF',N.float64,(4,1)),('C',N.float64,(4,4)),
                    ('Cesq',N.float64,(4,4)),('S',N.float64,(4,1)),
                    ('Sesq',N.float64,(4,1)),('NActive',N.float64,(4,1)),
                    ('activeEq',N.float64,(4,1)),('Nfree',N.intc),
                    ('freeS',numpy.intc))
                    
##typedef struct {
##  unsigned long sec[4] ;
##  double Y[4], Yesq[4], lambI, lambF ;
##  double C[4][4], Cesq[4][4] ;  #see 2darray.py for how to reference/dereference
##  double S[4], Sesq[4];
##  int Nactive, activeEq[4];
##  int Nfree, freeS[4];
##} PBdatapt;


class ConstraintEq(ctypes.Structure):
    _fields_=[("freeToConstrain",c_int),
            ("Nfree",c_int),
            ("freeS",c_int*4),
            ("Coeff",c_double*4)
             ]



##typedef struct {
##  int freeToConstrain; /* index 0-3 of free variable to constrain */
##  int Nfree, freeS[4];
##  double Coef[4];
##} ConstraintEq;


class He3CELL(ctypes.Structure):
    _fields_=[("name",c_char*64),
            ("tEmpty",c_double),
            ("tEmptySlope",c_double),
            ("L",c_double),
            ("D",c_double),
            ("R",c_double),
            ("nsL0",c_double),
            ("nsL0err",c_double),
            ("nsL",c_double),
            ("nsLerr",c_double)
             ]
#TODO check if there is a problem if a string is NULL terminated


##typedef struct {
##  char   name[64];
##  double tEmpty; /* at 1.77 A, may eventually have to put in wavelength dep */
##  double tEmptySlope; /* wavelength dependence */
##  double L;      /* max gas length in cm */
##  double D;      /* diameter in cm for id purposes only */
##  double R;      /* radius of curvature of windows, may want 2 values */
##  double nsL0;   /* uncorrected nsL at 1.77 A */
##  double nsL0err;
##  double nsL;    /* corrected for curvature */
##  double nsLerr;
##} He3CELL;       /* base He3 CELL constants */


class He3CELLpol(ctypes.Structure):
    _fields_=[("PHe",c_double),
            ("PHeErr",c_double),
            ("tEmptySlope",c_double),
            ("L",c_double),
            ("D",c_double),
            ("R",c_double),
            ("nsL0",c_double),
            ("nsL0err",c_double),
            ("nsL",c_double),
            ("nsLerr",c_double)
             ]

##typedef struct {
##  double PHe;    /* init He3 polarization */
##  double PHeErr; 
##  double hrsBeam; /* beamtime hours logged this exp */
##  double T;      /* cell decay time in hours in holding fld */
##  double Terr;
##  unsigned long startSecs; /* conv of startDate to seconds since Jan1 1971 */
##} He3CELLpol;

class expResol(ctypes.Structure):
    _fields_=[("waveRelWidth",c_double),
            ("angleVwidth",c_double),
            ("angleHwidth",c_double),
            ("double usedRadius",c_double)
             ]


##typedef struct {
##  double waveRelWidth;   /* std-dev-Lambda/Lambda */
##  double angleVwidth;    /* vertical angular resolution stnd dev */
##  double angleHwidth;    /* horizontal angular resolution stnd dev */
##  double usedRadius;     /* eff beam radius in cm for cell */
##} expResol;

class He3CELLexp(ctypes.Structure):
    _fields_=[("cell",He3CELL),
            ("pol",He3CELLpol),
            ("res",expResol)
             ]


##typedef struct {
##  He3CELL cell;
##  He3CELLpol pol;
##  expResol res;
##} He3CELLexp;

class efficiency(ctypes.Structure):
    _fields_=[("teff",c_double),
            ("terr",c_double),
            ("feff",c_double),
            ("ferr",c_double)
             ]



##typedef struct {
##  double teff;  /* transport efficiency */
##  double terr;
##  double feff;  /* flipper efficiency */
##  double ferr;
##} efficiency;

class PBsetup(ctypes.Structure):
    _fields_=[("P",He3CELLexp ),
            ("A",He3CELLexp ),
            ("eP",efficiency ),
            ("eA",efficiency )
             ]


##typedef struct {
##  He3CELLexp P;
##  He3CELLexp A;
##  efficiency eP;
##  efficiency eA;
##} PBsetup;



#pt=PBdatapt()
#for i in range(4):
#    pt.sec[i]=i

lib=ctypes.cdll['polarization.dll']


cell1=("PcellBT7Jan072008.txt")
cell2=("AcellBT7Jan72008.txt")
inputfile=("pol.out")
filename=("polarization.exe")
Argv=c_char_p*4

s=N.array([filename,cell1,cell2,inputfile])
argv=Argv(filename,cell1,cell2,inputfile)
#s=filename+"\0"+cell1+"\0"+cell2+"\0"+ inputfile
#mystr=ctypes.create_string_buffer(s)
#c_char_ap=ctypes.ARRAY(c_char_p,c_char)
#print repr(mystr.raw)
print type(s)
##mystrc = ctypes.cast(s, c_char_pp)
#data=N.ascontiguousarray(s)
#mystrc=data.ctypes.data_as(c_char_pp)
#print s
#print repr(mystrc.raw)
#print mystrc[:len(mystr)]
#mystr=ctypes.c_char_p(filename+" "+cell1+" "+cell2+" "+ inputfile)
#print mystr
#lib = N.ctypeslib.load_library('polarization.dll', '.')
polcor=lib.mypolcor

#polcor.argtyps=[c_int,c_char_pp]
polcor.restype=c_int
flag=polcor(4,ctypes.byref(argv))

#lib.PBcorrect(ctypes.byref(pt))
#print pt.sec[0]