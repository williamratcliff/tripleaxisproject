import copy
import sys

class Chair(object):
    def __init__(self,a=[],b=[]):
        print 'init'
        self._a=a
        self._b=b
        print 'a init',self._a,a
        print 'b init',self._b,b



if __name__=="__main__":
    spin=[0,1,0]
    myg=Chair()
    myg._a.append(copy.deepcopy(spin))
    print myg._a
    
    spin2=[1,0,0]
    bob=Chair()
    bob._a.append(spin2)
    print bob._a
    print sys.version