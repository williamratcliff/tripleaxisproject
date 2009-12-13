import numpy as N
import pylab
from matplotlib.ticker import NullFormatter, MultipleLocator
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator


def gen_as():
  xas=0.36
  aspos=N.array([[0.0000,   0.0000,1-xas],   
                 [0.0000,   0.0000,xas],    
                 [0.0000,   0.5000,.5-xas],    
                 [0.0000,   0.5000,.5+xas],       
                 [0.5000,   0.0000,.5-xas],  
                 [0.5000,   0.0000,.5+xas],    
                 [0.5000,   0.5000,1-xas],   
                 [0.5000,   0.5000,xas]],'Float64')
  return aspos



if __name__=="__main__":
  aspos=gen_as()
  xn=aspos[N.where(aspos[:,1]==0)[0],0]
  zn=aspos[N.where(aspos[:,1]==0)[0],2]
  xp=aspos[N.where(aspos[:,1]==0.5)[0],0]
  zp=aspos[N.where(aspos[:,1]==0.5)[0],2]
  
