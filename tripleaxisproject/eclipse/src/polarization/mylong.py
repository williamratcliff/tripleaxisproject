import numpy as N
import math
import ctypes
import cstruct
from ctypes import c_void_p, c_int, c_long, c_char, c_char_p,c_double,c_ulong
from ctypes import byref as _ref
c_void_pp = ctypes.POINTER(c_void_p)
c_char_pp = ctypes.POINTER(c_char_p)
c_int_p = ctypes.POINTER(c_int)
c_double_p=ctypes.POINTER(c_double)
c_ulong_p=ctypes.POINTER(c_ulong)
cellP="PcellBT7Jan072008.txt"
cellA="AcellBT7Jan72008.txt"



class indata(ctypes.Structure):
    _fields_=[("y",c_ulong_p)]

myin=indata()
dummy=N.array([10,23,37],'uint64')
myin.y=dummy.ctypes.data_as(c_ulong_p)
npts=dummy.shape[0]
mypolcorrect = N.ctypeslib.load_library('mylong.dll', '.')
mypolcorrect.PBcorrectData(ctypes.byref(myin),npts)
print dummy