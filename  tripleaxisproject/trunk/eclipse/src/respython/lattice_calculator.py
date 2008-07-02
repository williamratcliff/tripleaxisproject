import numpy as N
import pylab
import math
import unittest
from matplotlib.patches import Ellipse
eps=1e-3
pi=N.pi

def autovectorized(f):
     """Function decorator to do vectorization only as necessary.
     vectorized functions fail for scalar inputs."""
     def wrapper(input):
         if N.isscalar(input)==False:
             return N.vectorize(f)(input)
         return f(input)
     return wrapper



@autovectorized
def myradians(x):
    return math.radians(x)

#vecradians = N.vectorize(myradians, otypes=[double]) 

def sign(x):
    if x>0:
         ret=1
    if x<0:
        ret=-1
    if x==0:
        ret=0
    return ret

def blkdiag(g):
    "Returns a block diagonal matrix, given a list of square matrices, g"
    glen=len(g)
    n=0
    for i in range(glen):
        n=n+g[i].shape[0]


    gout=N.zeros((n,n))
    offset=0
    for i in range(glen):
        currblock=g[i]
        lenx,leny=currblock.shape
        for x in range(lenx):
            for y in range(leny):
                gout[x+offset,y+offset]=currblock[x,y]
        offset=offset+lenx
    return gout


def similarity_transform(A,B):
    G=N.dot(B,A.transpose())
    G2=N.dot(A,G)
    #
    return G2

class instrument:
    def __init__(self):
            self.tau_list={'pg(002)':1.87325, \
                            'pg(004)':3.74650, \
                            'ge(111)':1.92366, \
                            'ge(220)':3.14131, \
                            'ge(311)':3.68351, \
                            'be(002)':3.50702, \
                            'pg(110)':5.49806}
            return
    def get_tau(self,tau):
        return self.tau_list[tau]


class lattice:
    def __init__(self, a=N.array([2*N.pi, 2*N.pi], 'd'), \
                 b=N.array([2*N.pi, 2*N.pi], 'd'), \
                 c=N.array([2*N.pi, 2*N.pi], 'd'), \
                 alpha=N.array([90, 90], 'd'), \
                 beta=N.array([90, 90], 'd'), \
                 gamma=N.array([90, 90], 'd'), \
                 orient1=N.array([[1, 0, 0], [1, 0, 0]], 'd'), \
                 orient2=N.array([[0, 1, 0], [0, 1, 0]], 'd')):
        self.a=a
        self.b=b
        self.c=c
        self.alphad=alpha
        self.betad=beta
        self.gammad=gamma
        self.alpha=myradians(alpha)
        self.beta=myradians(beta)
        self.gamma=myradians(gamma)
        self.star()
        self.gtensor('lattice')
        self.gtensor('latticestar')
        self.npts=N.size(a)
        self.orient1=orient1.transpose()
        self.orient2=orient2.transpose()
        self.StandardSystem()
        self.instrument=instrument()
        return

    def scalar(self, x1, y1, z1, x2, y2, z2, lattice):
        "calculates scalar product of two vectors"
        if lattice=='lattice':
            a=self.a
            b=self.b
            c=self.c
            alpha=self.alpha
            beta=self.beta
            gamma=self.gamma
        if lattice=='latticestar':
            a=self.astar
            b=self.bstar
            c=self.cstar
            alpha=self.alphastar
            beta=self.betastar
            gamma=self.gammastar
        
        s=x1*x2*a**2+y1*y2*b**2+z1*z2*c**2+\
           (x1*y2+x2*y1)*a*b*N.cos(gamma)+\
           (x1*z2+x2*z1)*a*c*N.cos(beta)+\
           (z1*y2+z2*y1)*c*b*N.cos(alpha)
        return s

    def star(self):
        "Calculate unit cell volume, reciprocal cell volume, reciprocal lattice parameters"
        V=2*self.a*self.b*self.c*\
        N.sqrt(N.sin((self.alpha+self.beta+self.gamma)/2)*\
                N.sin((-self.alpha+self.beta+self.gamma)/2)*\
                N.sin((self.alpha-self.beta+self.gamma)/2)*\
                N.sin((self.alpha+self.beta-self.gamma)/2))
        self.Vstar=(2*N.pi)**3/V;
        self.astar=2*N.pi*self.b*self.c*N.sin(self.alpha)/V
        self.bstar=2*N.pi*self.a*self.c*N.sin(self.beta)/V
        self.cstar=2*N.pi*self.b*self.a*N.sin(self.gamma)/V
        self.alphastar=N.arccos((N.cos(self.beta)*N.cos(self.gamma)-\
                                 N.cos(self.alpha))/ \
                                (N.sin(self.beta)*N.sin(self.gamma)))
        self.betastar= N.arccos((N.cos(self.alpha)*N.cos(self.gamma)-\
                                 N.cos(self.beta))/ \
                                (N.sin(self.alpha)*N.sin(self.gamma)))
        self.gammastar=N.arccos((N.cos(self.alpha)*N.cos(self.beta)-\
                               N.cos(self.gamma))/ \
                              (N.sin(self.alpha)*N.sin(self.beta)))
        self.V=V
        return

    def angle2(self, x, y, z, h, k, l):
        "Calculate the angle between vectors in real and reciprocal space"
        "x,y,z are the fractional cell coordinates of the first vector,"
        "h,k,l are Miller indices of the second vector"
        phi=N.arccos(2*pi*(h*x+k*y+l*z)/self.modvec(x, y, z, 'lattice')/self.modvec(h, k, l, 'latticestar'))
        return phi


    def angle(self, x1, y1, z1, x2, y2, z2, lattice):
        "Calculate the angle between vectors in real and reciprocal space"
        "xi,yi,zi are the fractional cell coordinates of the vectors"
        phi=N.arccos(self.scalar(x1, y1, z1, x2, y2, z2, lattice)/self.modvec(x1, y1, z1, lattice)/self.modvec(x1, y1, z1, lattice))    
        return phi
    
    def modvec(self, x, y, z, lattice):
        "Calculates modulus of a vector defined by its fraction cell coordinates"
        "or Miller indexes"
        m=N.sqrt(self.scalar(x, y, z, x, y, z, lattice))
        return m

    def gtensor(self, lattice):
        "calculates the metric tensor of a lattice"
        g=N.zeros((3, 3, N.size(self.a)), 'd')
        #print 'shape ', g.shape
        if lattice=='lattice':
            a=self.a
            b=self.b
            c=self.c
            alpha=self.alpha
            beta=self.beta
            gamma=self.gamma
        if lattice=='latticestar':
            a=self.astar
            b=self.bstar
            c=self.cstar
            alpha=self.alphastar
            beta=self.betastar
            gamma=self.gammastar
        g[0, 0, :]=a**2;
        g[0, 1, :]=a*b*N.cos(gamma)
        g[0, 2, :]=a*c*N.cos(beta)

        g[1, 0, :]=g[0, 1, :]
        g[1, 1, :]=b**2
        g[1, 2, :]=c*b*N.cos(alpha)

        g[2, 0, :]=g[0, 2, :]
        g[2, 1, :]=g[1, 2, :]
        g[2, 2, :]=c**2
        if lattice=='lattice':
            self.g=g
        if lattice=='latticestar':
            self.gstar=g
        return

    def reciprocate(self, x, y, z, lattice):
        "calculate miller indexes of a vector defined by its fractional cell coords"
        if lattice=='lattice':
            g=self.g
        if lattice=='latticestar':
            g=self.gstar
        h=g[0, 0, :]*x+g[1, 0, :]*y+g[2, 0, :]*z;
        k=g[0, 1, :]*x+g[1, 1, :]*y+g[2, 1, :]*z;
        l=g[0, 2, :]*x+g[1, 2, :]*y+g[2, 2, :]*z;
        return h, k, l

    def vector(self, x1, y1, z1, x2, y2, z2, lattice):
        "calculates the fractional cell coordinates or Miller indexes of a vector"
        "product of two vectors, defined by their fractional cell coordinates or "
        "Miller idexes"
        if lattice=='lattice':
            g=self.gstar
            V=self.Vstar
        if lattice=='latticestar':
            g=self.g
            V=self.V
        g=g*V/(2*N.pi)**2
        x=y1*z2*g[0, 0, :]-z1*y2*g[0, 0, :]-x1*z2*g[1, 0, :]+z1*x2*g[1, 0, :]\
           +x1*y2*g[2, 0, :]-y1*x2*g[2, 0, :]
        y=y1*z2*g[0, 1, :]-z1*y2*g[0, 1, :]-x1*z2*g[1, 1, :]+z1*x2*g[1, 1, :]\
           +x1*y2*g[2, 1, :]-y1*x2*g[2, 1, :]
        z=y1*z2*g[0, 2, :]-z1*y2*g[0, 2, :]-x1*z2*g[1, 2, :]+z1*x2*g[1, 2, :]\
           +x1*y2*g[2, 2, :]-y1*x2*g[2, 2, :]
        return x,y,z

    def StandardSystem(self):
        orient1=self.orient1
        orient2=self.orient2
        modx=self.modvec(orient1[0, :], orient1[1, :], orient1[2, :], 'latticestar')
        x=N.copy(orient1)
        x[0, :]=x[0, :]/modx; # First unit basis vector
        x[1, :]=x[1, :]/modx;
        x[2, :]=x[2, :]/modx;
        
        proj=self.scalar(orient2[0, :], orient2[1, :], orient2[2, :], \
                    x[0, :], x[1, :], x[2, :], 'latticestar')
        
        y=N.copy(orient2)
        y[0, :]=y[0, :]-x[0, :]*proj; 
        y[1, :]=y[1, :]-x[1, :]*proj;
        y[2, :]=y[2, :]-x[2, :]*proj;

        mody=self.modvec(y[0, :], y[1, :], y[2, :], 'latticestar');

#    check for collinearity of orienting vectors

        try:
            if N.where(mody<=eps)[0].size>0:
                print 'ValueError'
                raise ValueError
            y[0, :]=y[0, :]/mody; # Second unit basis vector
            y[1, :]=y[1, :]/mody;
            y[2, :]=y[2, :]/mody;
    
            z=N.copy(y);
    
            z[0, :]=x[1, :]*y[2, :]-y[1, :]*x[2, :];
            z[1, :]=x[2, :]*y[0, :]-y[2, :]*x[0, :];
            z[2, :]=-x[1, :]*y[0, :]+y[1, :]*x[0, :];
    
            proj=self.scalar(z[0, :], z[1, :], z[2, :], x[0, :], x[1, :], x[2, :], 'latticestar');
    
            z[0, :]=z[0, :]-x[0, :]*proj; 
            z[1, :]=z[1, :]-x[1, :]*proj;
            z[2, :]=z[2, :]-x[2, :]*proj;
    
            proj=self.scalar(z[0, :], z[1, :], z[2, :], y[0, :], y[1, :], y[2, :], 'latticestar');
    
            z[0, :]=z[0, :]-y[0, :]*proj; 
            z[1, :]=z[1, :]-y[1, :]*proj;
            z[2, :]=z[2, :]-y[2, :]*proj;
    
            modz=self.modvec(z[0, :], z[1, :], z[2, :], 'latticestar');
    
            z[0, :]=z[0, :]/modz; #% Third unit basis vector
            z[1, :]=z[1, :]/modz;
            z[2, :]=z[2, :]/modz;     
            
            self.x=x
            self.y=y
            self.z=z    
        except ValueError:
            print 'ORIENTATION VECTORS ARE COLLINEAR x,y,z not set'    
        return
    
    def S2R(self, qx, qy, qz):
        "Given cartesian coordinates of a vector in the S System, calculate its Miller indexes."
        x=self.x
        y=self.y
        z=self.z
        H=qx*x[0, :]+qy*y[0, :]+qz*z[0, :];
        K=qx*x[1, :]+qy*y[1, :]+qz*z[1, :];
        L=qx*x[2, :]+qy*y[2, :]+qz*z[2, :];
        q=N.sqrt(qx**2+qy**2+qz**2);
        return H, K, L, q
    
    def R2S(self, H, K, L):
        "Given reciprocal-space coordinates of a vecotre, calculate its coordinates in the Cartesian space."
        x=self.x
        y=self.y
        z=self.z
        qx=self.scalar(H, K, L, x[0, :], x[1, :], x[2, :], 'latticestar');
        qy=self.scalar(H, K, L, y[0, :], y[1, :], y[2, :], 'latticestar');
        qz=self.scalar(H, K, L, z[0, :], z[1, :], z[2, :], 'latticestar');
        q=self.modvec(H, K, L, 'latticestar');
        return qx, qy, qz, q

    def SpecWhere(self,M2,S1,S2,A2,EXP):
        """ For given values of M3,S1,S2 and A2 spectrometer motors (AKA M2,M3,M4 and M6)
        and spectrometer and sample parameters specified in EXP calculates the wave vector
        transfer in the sample (H, K, L), Q=|(H,K,L)|, energy tranfer E, and incident
        and final neutron energies."""


        taum=N.empty(M2.shape,'d')
        taua=N.empty(M2.shape,'d')
        for ind in range(M2.shape[0]):
             taum[ind]=self.instrument.get_tau(EXP[ind]['mono']['tau'])
        for ind in range(M2.shape[0]):
             taua[ind]=self.instrument.get_tau(EXP[ind]['ana']['tau'])
        ki=taum/N.sqrt(2.0-2*N.cos(M2))
        Ei=2.072142*ki**2
        kf=taua/N.sqrt(2.0-2*N.cos(A2))
        Ef=2.072142*kf**2
        E=Ei-Ef
        Q=N.sqrt(ki**2+kf**2-2*ki*kf*N.cos(S2))
        orienta=self.x
        orientb=self.y
        #phi=-atan2(-kf.*sin(S2), ki-kf.*cos(S2)); %Angle from ki to Q
        delta=N.absolute(N.arccos( (Q**2+ki**2-kf**2)/(2*ki*Q)))
        psi=S1+delta-pi/2 #Angle from first orienting vector to to Q
        qx=Q*N.cos(psi)
        qy=Q*N.sin(psi)
        H=qx*orienta[0]+qy*orientb[0]
        K=qx*orienta[1]+qy*orientb[1]
        L=qx*orienta[2]+qy*orientb[2]
        return H,K,L,E,Q,Ei,Ef

       
class TestLattice(unittest.TestCase):

    def setUp(self):
        a=N.array([2*pi],'d')
        b=N.array([2*pi],'d')
        c=N.array([2*pi],'d')
        alpha=N.array([90],'d')
        beta=N.array([90],'d')
        gamma=N.array([90],'d')
        orient1=N.array([[1,0,0]],'d')
        orient2=N.array([[0,1,1]],'d')
        self.fixture = lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                               orient1=orient1,orient2=orient2)
    
    def test_astar(self):
        self.assertAlmostEqual(self.fixture.astar[0],1.0,2,'astar Not equal to '+str(1.0))
    def test_bstar(self):
        self.assertAlmostEqual(self.fixture.bstar[0],1.0,2,'bstar Not equal to '+str(1.0))
    def test_cstar(self):
        self.assertAlmostEqual(self.fixture.cstar[0],1.0,2,'cstar '+str(self.fixture.cstar[0])+' Not equal to '+str(1.0))
    def test_alphastar(self):
        self.assertAlmostEqual(self.fixture.alphastar[0],pi/2,2,'alphastar Not equal to '+str(pi/2))
    def test_betastar(self):
        self.assertAlmostEqual(self.fixture.betastar[0],pi/2,2,'betastar Not equal to '+str(pi/2))
    def test_gammastar(self):
        self.assertAlmostEqual(self.fixture.gammastar[0],pi/2,2,'gammastar Not equal to '+str(pi/2))
    def test_V(self):
        self.assertAlmostEqual(self.fixture.V[0],248.0502,2,'V Not equal to '+str(248.0502))
    def test_Vstar(self):
        self.assertAlmostEqual(self.fixture.Vstar[0],1.0,2,'Vstar Not equal to '+str(1.0))
    def test_g(self):
        #print self.fixture.g
        self.assertAlmostEqual((self.fixture.g[:,:,0][0,0]),39.4784*(N.eye(3)[0,0]) ,2,'g Not equal to '+str(39.4784 ))
    def test_gstar(self):
        #print self.fixture.gstar
        self.assertAlmostEqual(self.fixture.gstar[:,:,0][0,0],1.0*N.eye(3)[0,0] ,2,'gstar Not equal to '+str(1.0 ))
 
    def test_StandardSystem_x(self):
 #       #print self.fixture.gstar
        self.assertAlmostEqual(self.fixture.x[0],1.0 ,2,'Standard System x Not equal to '+str(1.0 ))
 
    
      
                               
#    def test_zeroes(self):
#        self.assertEqual(0 + 0, 0)
#        self.assertEqual(5 + 0, 5)
#        self.assertEqual(0 + 13.2, 13.2)
#
#    def test_positive(self):
#        self.assertEqual(123 + 456, 579)
#        self.assertEqual(1.2e20 + 3.4e20, 3.5e20)
#
#    def test_mixed(self):
#        self.assertEqual(-19 + 20, 1)
#        self.assertEqual(999 + -1, 998)
#        self.assertEqual(-300.1 + -400.2, -700.3)
#        

if __name__=="__main__":
    if 1:
        a=N.array([2*pi],'d')
        b=N.array([8],'d')
        c=N.array([11],'d')
        alpha=N.array([87],'d')
        beta=N.array([52],'d')
        gamma=N.array([100],'d')
        orient1=N.array([[0,1,0]],'d')
        orient2=N.array([[1,0,0]],'d')
        mylattice=lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                               orient1=orient1,orient2=orient2)
        H=N.array([1],'d');K=N.array([0],'d');L=N.array([0],'d');W=N.array([0],'d')
        EXP={}
        EXP['ana']={}
        EXP['ana']['tau']='pg(002)'
        EXP['mono']={}
        EXP['mono']['tau']='pg(002)';
        EXP['ana']['mosaic']=30
        EXP['mono']['mosaic']=30
        EXP['sample']={}
        EXP['sample']['mosaic']=10
        EXP['sample']['vmosaic']=10
        EXP['hcol']=N.array([40, 10, 20, 80],'d')
        EXP['vcol']=N.array([120, 120, 120, 120],'d')
        EXP['infix']=-1 #positive for fixed incident energy
        EXP['efixed']=14.7
        EXP['method']=0
        setup=[EXP]  
        M2=myradians(N.array([41.177]))
        A2=myradians(N.array([41.177]))
        S1=myradians(N.array([66.4363]))
        S2=myradians(N.array([37.6547]))
        H,K,L,E,Q,Ei,Ef=mylattice.SpecWhere(M2,S1,S2,A2,setup)
        print 'H ',H
        print 'K ',K
        print 'L ',L
        print 'E ',E
        print 'Q ',Q
        print 'Ei ',Ei
        print 'Ef ',Ef
##        print 'M2 ',M2
##        print 'A2 ',A2
##        print 'S1 ',S1
##        print 'S2 ',S2

#    unittest.main()