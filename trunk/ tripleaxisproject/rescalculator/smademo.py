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
    Deltax=p[0] #					% Gap at the AF zone-center in meV for x-axis mode
    Deltay=p[1]	#				% Gap at the AF zone-center in meV for y-axis mode
    Deltaz=p[2]	#				% Gap at the AF zone-center in meV for z-axis mode
    cc=p[3] #						% Bandwidth in meV
    Gamma=p[4] #					    % Intrinsic HWHM of exccitation in meV
    #I=p[6);							% Intensity prefactor It will be used in the PrefDemo function, not here!
    #bgr=p[7);						% Background. It will be used in the PrefDemo function, not here!

    # Calculate the dispersion
    omegax=sqrt(cc**2*(sin(2*pi*H))**2+Deltax**2)
    omegay=sqrt(cc**2*(sin(2*pi*H))**2+Deltay**2)
    omegaz=sqrt(cc**2*(sin(2*pi*H))**2+Deltaz**2)
    w0[0,:]=omegax
    w0[1,:]=omegay
    w0[2,:]=omegaz

    # Intensity scales as (1-cos(2*pi*H))/omega0 for each of the three modes:
    sqw=zeros((3,size(H)),'float64')
    #print 'sqw ',sqw.shape
    #print 'lorx ',lorx.shape
    #print 'done'
    sqw[0,:]=lorx*(1-cos(pi*H))/omegax/2
    sqw[1,:]=lory*(1-cos(pi*H))/omegay/2
    sqw[2,:]=lorz*(1-cos(pi*H))/omegaz/2
    #print 'sqw ',sqw
    #% Now set all energy widths of all branches to Gamma
    HWHM=ones(sqw.shape,'float64')*Gamma
    return w0,sqw,HWHM
