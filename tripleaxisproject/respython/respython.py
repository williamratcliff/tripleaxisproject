import numpy as N
import math
import pylab

class instrument:
    def __init__(self,tau):
            self.tau_list={'pg(002)':1.87325, \
                            'pg(004)':3.74650,\
                            'ge(111)':1.92366,\
                            'ge(220)':3.14131,\
                            'ge(311)':3.68351,\
                            'be(002)':3.50702,\
                            'pg(110)':5.49806}
            self.tau=self.tau_list[tau]
            return
    def get_tau(self):
        return self.tau

class lattice:
    def __init__(self,a=2*N.pi,b=2*N.pi,c=2*N.pi,alpha=90,beta=90,gamma=90):
        self.a=a
        self.b=b
        self.c=c
        self.alphad=alpha
        self.betad=beta
        self.gammad=gamma
        self.alpha=math.radians(alpha)
        self.beta=math.radians(beta)
        self.gamma=math.radians(gamma)
        self.star()
        self.gtensor('lattice')
        self.gtensor('latticestar')
        return

    def scalar(self,x1,y1,z1,x2,y2,z2,lattice):
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

    def angle2(self,x,y,z,h,k,l):
        "Calculate the angle between vectors in real and reciprocal space"
        "x,y,z are the fractional cell coordinates of the first vector,"
        "h,k,l are Miller indices of the second vector"
        phi=N.arccos(2*pi*(h*x+k*y+l*z)/self.modvec(x,y,z,'lattice')/modvec(h,k,l,'latticestar'))
        return phi


    def angle(self,x1,y1,z1,x2,y2,z2,lattice):
        "Calculate the angle between vectors in real and reciprocal space"
        "xi,yi,zi are the fractional cell coordinates of the vectors"
        phi=N.arccos(self.scalar(x1,y1,z1,x2,y2,z2,lattice)/self.modvec(x1,y1,z1,lattice)/self.modvec(x1,y1,z1,lattice))    
        return phi
    
    def modvec(self,x,y,z,lattice):
        "Calculates modulus of a vector defined by its fraction cell coordinates"
        "or Miller indexes"
        m=N.sqrt(self.scalar(x,y,z,x,y,z,lattice))
        return m

    def gtensor(self,lattice):
        "calculates the metric tensor of a lattice"
        g=N.zeros((3,3,N.size(self.a)),'d')
        print 'shape ', g.shape
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
        g[0,0,:]=a**2;
        g[0,1,:]=a*b*N.cos(gamma)
        g[0,2,:]=a*c*N.cos(beta)

        g[1,0,:]=g[0,1,:]
        g[1,1,:]=b**2
        g[1,2,:]=c*b*N.cos(alpha)

        g[2,0,:]=g[0,2,:]
        g[2,1,:]=g[1,2,:]
        g[2,2,:]=c**2
        if lattice=='lattice':
            self.g=g
#            print 'lattice'
        if lattice=='latticestar':
            self.gstar=g
#            print 'latticestar'
#        print g
        return

    def reciprocate(self,x,y,z,lattice):
        "calculate miller indexes of a vector defined by its fractional cell coords"
        if lattice=='lattice':
            g=self.g
        if lattice=='latticestar':
            g=self.gstar
        h=g(1,1,:)*x+g(2,1,:)*y+g(3,1,:)*z;
        k=g(1,2,:)*x+g(2,2,:)*y+g(3,2,:)*z;
        l=g(1,3,:)*x+g(2,3,:)*y+g(3,3,:)*z;
        return h,k,l

    def vector(self,x1,y1,z1,x2,y2,z2,lattice):
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
        x=y1*z2*g(1,1,:)-z1*y2*g(1,1,:)-x1*z2*g(2,1,:)+z1*x2*g(2,1,:)\
           +x1*y2*g(3,1,:)-y1*x2*g(3,1,:)
        y=y1*z2*g(1,2,:)-z1*y2*g(1,2,:)-x1*z2*g(2,2,:)+z1*x2*g(2,2,:)\
           +x1*y2*g(3,2,:)-y1*x2*g(3,2,:)
        z=y1*z2*g(1,3,:)-z1*y2*g(1,3,:)-x1*z2*g(2,3,:)+z1*x2*g(2,3,:)\
           +x1*y2*g(3,3,:)-y1*x2*g(3,3,:)
 



if __name__=="__main__":
    mylattice=lattice();
    x1=N.array([1.0,1.0],'d'); y1=N.array([1.0,1.0],'d'); z1=N.array([1.0,1.0],'d'); x2=x1; y2=y1; z2=z1;
#    print 'scalar ', mylattice.scalar(x1,y1,z1,x2,y2,z2,'lattice')
#    print 'me ',mylattice.gtensor('lattice')
    myinstrument=instrument('pg(004)');
    print myinstrument.get_tau()
