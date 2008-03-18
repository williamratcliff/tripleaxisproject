from numpy import *

def SMADemo(H,K,L,p):
    """This is an example of a cross section function for use with ConvRes.m.
    This particular function calculates the cross section for a gapped
    excitations in a 1-dimensional antiferromagnet. "a" is the chain axis.
    The polarization factors for each mode are NOT calculated here, but
    should be included in the prefactor function instead. This function is
    meant to be used together with the prefactor function PrefDemo.m
    Arguments H-W and are vectors, so dont forget to use ".*" instead of "*", etc."""


    #% Extract the three parameters contained in "p":
    hcenter=p[0]
    kcenter=p[1]
    lcenter=p[2]
    corrh=p[3]
    corrk=p[3]
    corrl=p[3]

    # Calculate the dispersion
    w0=zeros((1,size(H)),'float64')

    # Intensity scales as (1-cos(2*pi*H))/omega0 for each of the three modes:
    #sqw=zeros((3,size(H)),'float64')
    #print 'sqw ',sqw.shape
    #print 'lorx ',lorx.shape
    #print 'done'
    qx=H-hcenter
    qy=K-kcenter
    qz=L-lcenter
    num=1.0*corrh*corrk*corrl/pi**2
    den=(qx**2*corrh**2+qy**2*corrk**2+qz**2*corrl**2+1)
    L=num/den**2.0
    sqw=L
    #print 'sqw ',sqw
    #% Now set all energy widths of all branches to Gamma
    HWHM=zeros(w0.shape,'float64')
    return w0,sqw,HWHM
