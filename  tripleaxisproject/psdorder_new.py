import numpy as N
import pylab
import scriptutil as SU
import re
import readncnr2 as readncnr
import simple_combine
#import scipy
from scipy.optimize import leastsq
import copy
import scipy.odr
pi=N.pi


if __name__=='__main__':
        mydirectory=r'C:\13165\13165\data'
        myfilebase='PSD_A4_SpacingApr1108.dat'
        myend='bt7'
        flist=[]
        rlist=[56421,56420,56419]
        myfilestr=r'C:\13165\PSD_A4_SpacingApr1108.dat'
        a4_pos=N.loadtxt(myfilestr, unpack=True)
        print a4_pos.T.flatten()
        myfilestr=r'C:\13165\PSD_A4_SpacingApr1108.dat'
        ch_eff=N.loadtxt(myfilestr, unpack=True)
        print a4_pos.T.flatten()

