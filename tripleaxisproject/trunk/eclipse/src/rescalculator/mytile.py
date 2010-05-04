import numpy as N

#def mytile(h,k):
#    """takes input arrays h,k and returns a grid""


if __name__=='__main__':
    ht=N.linspace(-1,1,21)
    kt=N.linspace(-2,2,41)
    h,k=N.meshgrid(ht,kt)
    print h.flatten(),k.flatten()