import numpy as N
import dct


if __name__=="__main__":
    
    myfilestr=r'c:\structfactors.dat'
    h,k,l,chi,sig=N.loadtxt(myfilestr)
    