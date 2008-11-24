import numpy as N
import sympy
from sympy import cos,exp,I,sin,latex,pngview


class atom:
    def __init__(self,spin=[0,0,1],pos=[0,0,0],neighbors=[],interactions=[],label=0,Dx=0,Dy=0,Dz=0,cell=0,int_cell=[]):
        self.spin=spin
        self.pos=N.array(pos,'float64')
        self.neighbors=neighbors
        self.interactions=interactions
        self.label=label
        self.Dx=Dx
        self.Dy=Dy
        self.Dz=Dz
        self.cell=cell
        self.int_cell=[]


def generate_sabn(N_atoms):
    "generate spins in local coordinate system, with Z as quantization axis"
    Sabn=[]
    S=sympy.Symbol("S",real=True)
    for i in range(N_atoms):
        c=sympy.Symbol('c%d'%(i,),commutative=False)
        cd=sympy.Symbol('cd%d'%(i,),commutative=False)
        curr=sympy.matrices.Matrix([sympy.sqrt(S/2.0)*(c+cd),sympy.sqrt(S/2.0)*(c-cd)/I,S-cd*c])
        Sabn.append(curr.reshape(3,1))
    return Sabn


def generate_sxyz(Sabn,atomlist):
    "transform spins from local coordinate system to global system"
    Sxyz=[]
    i=0
    for currS in Sabn:
        tempS=atomlist[i].spin*currS
        tempS=tempS.reshape(1,3)
        Sxyz.append(tempS)
        i=i+1
    return Sxyz


def generate_atoms_rot():
    "generate some atoms for our test case"
    if 1:
        "make 1st atom"
        spin0=sympy.matrices.Matrix([[1,0,0],[0,1,0],[0,0,1]])
        pos0=[0,0,0]
        neighbors=[1]
        interactions=[0]
        atom0=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions)
        "make 2nd atom"
        pos0=[1,0,0]
        P=sympy.Symbol('P',real=True,commutative=True)
        
        spin0=sympy.matrices.Matrix([[cos(P),-sin(P),0],[sin(P),cos(P),0],[0,0,1]])
        print 'spin0 converted',spin0
        neighbors=[0]
        interactions=[0]
        atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions)
        atomlist=[atom0,atom1]       
    return atomlist


def generate_hdef(atom_list,Jij,Sxyz,N_atoms_uc,N_atoms):
    N_atoms=len(atom_list)
    Hdef=0
    print 'Jij',len(Jij)
    
    
    #sympy.matrices.Matrix.multiply
    #for i in range(N_atoms):
    for i in range(N_atoms_uc): #correct
        N_int=len(atom_list[i].interactions)
        for j in range(N_int):
        #            currS_transpose=N.reshape(currS,(3,1))
        #tempS=atomlist[i].spin*currS_transpose
        #tempS=N.array(tempS)
        #tempS=N.ravel(tempS)
            
            print 'i',i,'j',j
            if 0:
                Hij=N.matrix(Sxyz[i])*atom_list[i].spin.T
                print 'making Ham'
                print 'spin i', atom_list[i].spin.T
                print 'Sxyz i', N.matrix(Sxyz[i]),N.matrix(Sxyz[i]).shape
                print 'Hijtemp',Hij
                Hij=Hij*Jij[atom_list[i].interactions[j]]
                print 'Jij', Jij[atom_list[i].interactions[j]]
                print 'Hij*Jij', Hij
                Hij=Hij*atom_list[atom_list[i].neighbors[j]].spin#
                print 'Hijtemp3', Hij.shape
            if 1:
                print 'Sxyz i', N.matrix(Sxyz[i]),N.matrix(Sxyz[i]).shape
                Hij=Sxyz[i]*Jij[atom_list[i].interactions[j]]
                print 'S*Jij', Hij,Hij.shape
            Sxyz_transpose=Sxyz[atom_list[i].neighbors[j]].reshape(3,1)
            #Sxyz_transpose=Sxyz_transpose.(3,1))
            print 'Sxyz.T',Sxyz_transpose.shape
            print 'Hij before multiply', Hij, Hij.shape
            Hij=Hij*Sxyz_transpose
            print 'Hij*Sxyz.T',Hij,Hij.shape
            Hij=Hij[0]
            #Hij=Hij+myterm
            Hij=-Hij-atom_list[i].Dx*Sxyz[i][0]**2-atom_list[i].Dy*Sxyz[i][1]**2-atom_list[i].Dz*Sxyz[i][2]**2
            Hdef=Hdef+Hij
    print 'generated hdef'
    print Hdef,Hdef.atoms(sympy.Symbol)
    return Hdef




def fouriertransform(atom_list,Jij,Hlin,k,N_atoms_uc,N_atoms):
    #N_atoms=len(atom_list)
    #N_atoms_uc=1
    #N_atoms_uc=N_atoms
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    print 'fourier'
    print Hlin
    print Hlin.atoms(sympy.Symbol)
    print 'expand'
    #Hlin=Hlin.expand()
    print Hlin.atoms(sympy.Symbol)
    #print Hlin
    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol('c%d'%(i,),commutative=False,real=True)
        cdi=sympy.Symbol('cd%d'%(i,),commutative=False,real=True)
        cki=sympy.Symbol('ck%d'%(i,),commutative=False,real=True)
        ckdi=sympy.Symbol('ckd%d'%(i,),commutative=False,real=True)
        cmki=sympy.Symbol('cmk%d'%(i,),commutative=False,real=True)
        cmkdi=sympy.Symbol('cmkd%d'%(i,),commutative=False,real=True)
        ri=atom_list[i].pos
        N_int=len(atom_list[i].interactions)
        for j in range(N_atoms):
            rj=atom_list[j].pos
            j2=i#atom_list[i].neighbors[j]
            cj=sympy.Symbol('c%d'%(j,),commutative=False,real=True)
            cdj=sympy.Symbol('cd%d'%(j,),commutative=False,real=True)
            ckj=sympy.Symbol('ck%d'%(j2,),commutative=False,real=True)
            ckdj=sympy.Symbol('ckd%d'%(j2,),commutative=False,real=True)
            cmkj=sympy.Symbol('cmk%d'%(j2,),commutative=False,real=True)
            cmkdj=sympy.Symbol('cmkd%d'%(j2,),commutative=False,real=True)
            diffr=ri-rj
            kmult=N.dot(k,diffr)
            t1=1.0/2*(ckdi*cmkdj*exp(-I*kmult)+cmkdi*ckdj*exp(I*kmult)               )
            t2=1.0/2*(cki*cmkj*exp(I*kmult)+cmki*ckj*exp(-I*kmult))
            t3=1.0/2*(ckdi*ckj*exp(-I*kmult)+cmkdi*cmkj*exp(I*kmult))
            t4=1.0/2*(cki*ckdj*exp(I*kmult)+cmki*cmkdj*exp(-I*kmult))
            t5=1.0/2*(ckdj*ckj+cmkdj*cmkj)
            print 'i',i,'j',j
            print 'ci',ci,'cj',cj,'cdi',cdi,'cdj',cdj
            print 't1',t1
            print 't2',t2
            print 't3',t3
            print 't4',t4
            print 't5',t5
            a1=sympy.Symbol('a1')
            a2=sympy.Symbol('a2')
            a3=sympy.Symbol('a3')
            a4=sympy.Symbol('a4')
            a5=sympy.Symbol('a5')
            f1=cdi*cdj
            print 'f1',f1,f1.atoms(sympy.Symbol)
            Hlin=Hlin.subs(f1,t1)
            print 'H1',Hlin,Hlin.atoms(sympy.Symbol)
            f2=ci*cj
            print 'f2',f2,f2.atoms(sympy.Symbol)
            Hlin=Hlin.subs(f2,t2)
            print 'H2',Hlin,Hlin.atoms(sympy.Symbol)
            f3=cdi*cj
            print 'f3',f3,f3.atoms(sympy.Symbol)
            Hlin=Hlin.subs(f3,t3)
            print 'H3',Hlin,Hlin.atoms(sympy.Symbol)
            f4=ci*cdj
            print 'f4',f4,f4.atoms(sympy.Symbol)
            Hlin=Hlin.subs(f4,t4)
            print 'H4',Hlin,Hlin.atoms(sympy.Symbol)
            f5=cdj*cj
            print 'f5',f5,f5.atoms(sympy.Symbol)
            Hlin=Hlin.subs(f5,t5)
            print 'H5',Hlin,Hlin.atoms(sympy.Symbol)
    #print t1
    return Hlin#.expand()        
        
        
if __name__=="__main__":
    atom_list=generate_atoms_rot()
    N_atoms=2
    Sabn=generate_sabn(N_atoms)        
    Sxyz=generate_sxyz(Sabn,atom_list)
    N_atoms_uc=1
    Jij=[sympy.matrices.Matrix([[1,0,0],[0,1,0],[0,0,1]])]
    Hdef=generate_hdef(atom_list,Jij,Sxyz,N_atoms_uc,N_atoms)
    print 'Hdef',Hdef
    print Hdef.atoms()
    
    #pngview(Hdef)
    #print latex(Hdef)
    #print_matplotlib(latex(Hdef)) 
    kx=sympy.Symbol('kx',real=True)
    ky=sympy.Symbol('ky',real=True)
    kz=sympy.Symbol('kz',real=True)
    k=[kx,ky,kz]
    Hlin=Hdef
    Hfou=fouriertransform(atom_list,Jij,Hlin,k,N_atoms_uc,N_atoms)
    