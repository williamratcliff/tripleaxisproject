import sympy
import numpy as N


class atom:
    def __init__(self,spin=[0,0,1],pos=[0,0,0],neighbors=[],interactions=[],label=0,Dx=0,Dy=0,Dz=0):
        self.spin=spin
        self.pos=pos
        self.neighbors=neighbors
        self.interactions=interactions
        self.label=label
        self.Dx=Dx
        self.Dy=Dy
        self.Dz=Dz

def generate_sabn(N_atoms):
    Sabn=[]
    S=sympy.Symbol("S")
    for i in range(N_atoms):
        c=sympy.Symbol("c%d"%(i,),commutative=False)
        cd=sympy.Symbol("cd%d"%(i,),commutative=False)
        curr=[sympy.sqrt(S/2)*(c+cd),sympy.sqrt(S/2)*(c-cd)/sympy.I,S-cd*c]
        Sabn.append(curr)
    return Sabn

def generate_sxyz(N_atoms):
    Sxyz=[]
    for i in range(N_atoms):
        S=sympy.Symbol("Sxyz%d"%(i,),commutative=False)
        Sxyz.append(S)
    return Sxyz

def generate_hdef(atom_list,Jij,Sxyz):
    N_atoms=len(atom_list)
    Hdef=0
    print 'atom_list',atom_list
    print 'Sxyz',Sxyz
    print 'Jij',Jij
    for i in range(N_atoms):
        N_int=len(atom_list[i].interactions)
        for j in range(N_int):
            #print Jij[atom_list[i].interactions[j]]
            #print '1',Sxyz[i]
            #print Sxyz[atom_list[i].neighbors[j]]
            Hij=Sxyz[i]*Jij[atom_list[i].interactions[j]]#
            #print type(Hij)
            Sxyz_transpose=N.matrix(Sxyz[atom_list[i].neighbors[j]])
            Sxyz_transpose=N.reshape(Sxyz_transpose,(3,1))
            Hij=Hij*Sxyz_transpose
            Hij=Hij-atom_list[i].Dx*Sxyz[i][0]**2-atom_list[i].Dy*Sxyz[i][1]**2-atom_list[i].Dz*Sxyz[i][2]**2
            Hdef=Hdef+Hij
    return Hdef[0][0]

def generate_atoms():
    spin0=[0,0,1]
    pos0=[0,0,0]
    neighbors=[1]
    interactions=[0]
    atom0=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=0)
    pos0=[0,0.5,0]
    neighbors=[0]
    interactions=[0]
    atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1)
    atomlist=[atom0,atom1]
    return atomlist
   
    
if __name__=='__main__':

    if 0:
        A = sympy.Symbol("A", commutative=False)
        B = sympy.Symbol("B", commutative=False)
        e = (A+B)**2
        g=e.expand()
        l=g.subs('A*B',1)
        print l
    if 0:
        x = sympy.Symbol('x')
        p = sympy.Wild('p',exclude='x')
        q = sympy.Wild('q',exclude='x')
        myexp=(5*x**2 + 3*x)
        r=myexp.match(p*x**2+q*x )
        print r
        print 1./20
    if 1:
        N_atoms=3
        Sabn=generate_sabn(N_atoms)
        print Sabn[0]
        Sxyz=generate_sxyz(N_atoms)
        atom_list=generate_atoms()
        Jij=[N.matrix([[1,0,0],[0,1,0],[0,0,1]])]
        Hdef=generate_hdef(atom_list,Jij,Sabn)
        S = sympy.Symbol('S')
        Hdef=Hdef.expand()
        Hdef=Hdef.as_poly(S)
        p = sympy.Wild('p',exclude='S')
        q = sympy.Wild('q',exclude='S')
        r = sympy.Wild('r',exclude='S')
        l = sympy.Wild('l',exclude='S')
        #Hlin=Hdef.match(p*S**2+S*q+l)
        print Hdef
        print Hdef.coeffs[0]
        print Hdef.coeffs[1]
        Hlin=Hdef.coeffs[0]*S**2+Hdef.coeffs[1]*S
        print Hlin