import sympy
import numpy as N
#I=sympy.I
#I=1.0j
#pi=sympy.pi
from sympy import exp,I,pi,sin,cos
import pylab
#import readfiles
from sympy import pngview,latex
import scipy.linalg


def print_matplotlib(s):
    pylab.text(0,0,s)
    pylab.axis('off')
    pylab.show()
    return 