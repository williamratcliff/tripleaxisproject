import numpy as N
import sympy
from sympy import cos,exp,I,sin,latex,pngview


class atom:
    def __init__(self,spin=[0,0,1],pos=[0,0,0],neighbors=[],interactions=[],label=0,Dx=0,Dy=0,Dz=0,cell=0,int_cell=[]):
        self.spin=spin
        self.pos=sympy.matrices.Matrix(pos)
        self.neighbors=neighbors
        self.interactions=interactions
        self.label=label
        self.Dx=Dx
        self.Dy=Dy
        self.Dz=Dz
        self.cell=cell
        self.int_cell=[]




def generate_hdef(N_atoms_uc):
    Jij=[sympy.matrices.Matrix([[1,0,0],[0,1,0],[0,0,1]])]
    S=sympy.Symbol("S")
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
        W=sympy.Symbol('W',real=True,commutative=True)
        X=sympy.Symbol('X',real=True,commutative=True)
        spin0=sympy.matrices.Matrix([[W,-X,0],[W,X,0],[0,0,1]])
        #spin0=sympy.matrices.Matrix([[cos(P),-sin(P),0],[sin(P),cos(P),0],[0,0,1]])
        spin0=sympy.matrices.Matrix([[1,0,0],[0,1,0],[0,0,1]])
        neighbors=[0]
        interactions=[0]
        atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions)
        atom_list=[atom0,atom1]       

    
    
    N_atoms=len(atom_list)
    Hdef=0
    Sabn=[]
    
    if 0:
        for i in range(N_atoms):
            ci=sympy.Symbol('c%d'%(i,),commutative=False)
            cdi=sympy.Symbol('cd%d'%(i,),commutative=False)
            curr=sympy.matrices.Matrix([sympy.sqrt(S/2)*(ci+cdi),sympy.sqrt(S/2)*(ci-cdi)/I,S-cdi*ci])
            print 'Sabn',curr
            Sabn.append(curr.reshape(3,1))
           
    c0=sympy.Symbol('c0',commutative=False,real=True)
    cd0=sympy.Symbol('cd0',commutative=False,real=True)
    c1=sympy.Symbol('c1',commutative=False,real=True)
    cd1=sympy.Symbol('cd1',commutative=False,real=True)
    curr=sympy.matrices.Matrix([sympy.sqrt(S/2)*(c0+cd0),sympy.sqrt(S/2)*(c0-cd0)/I,S-cd0*c0])
    Sabn.append(curr.reshape(3,1))
    curr=sympy.matrices.Matrix([sympy.sqrt(S/2)*(c1+cd1),sympy.sqrt(S/2)*(c1-cd1)/I,S-cd1*c1])
    Sabn.append(curr.reshape(3,1))
    Sxyz=[]
    i=0
    for currS in Sabn:
        tempS=atom_list[i].spin*currS
        #tempS=tempS.reshape(1,3)
        Sxyz.append(tempS.T)
        i=i+1
    for i in range(N_atoms_uc): #correct
        N_int=len(atom_list[i].interactions)
        for j in range(N_int):            
            print 'i',i,'j',j
            if 1:
                print 'Sxyz i', Sxyz[i],Sxyz[i].shape
                Hij=sympy.matrices.Matrix.multiply(Sxyz[i],Jij[atom_list[i].interactions[j]])
                print 'S*Jij', Hij,Hij.shape
            Sxyz_transpose=Sxyz[atom_list[i].neighbors[j]].T#reshape(3,1)
            #Sxyz_transpose=Sxyz_transpose.(3,1))
            print 'Sxyz.T',Sxyz_transpose.shape
            print 'Hij before multiply', Hij, Hij.shape
            #Hij=Hij*Sxyz_transpose
            Hij=sympy.matrices.Matrix.multiply(Hij,Sxyz_transpose)
            print 'Hij*Sxyz.T',Hij,Hij.shape
            Hij=Hij[0,0]
            Hij=-Hij-atom_list[i].Dx*Sxyz[i][0]**2-atom_list[i].Dy*Sxyz[i][1]**2-atom_list[i].Dz*Sxyz[i][2]**2
            Hdef=Hdef+Hij
    print 'generated hdef'
    print Hdef,Hdef.atoms(sympy.Symbol)

    Hlin=Hdef
 
    kx=sympy.Symbol('kx',real=True)
    ky=sympy.Symbol('ky',real=True)
    kz=sympy.Symbol('kz',real=True)
    k=sympy.matrices.Matrix([kx,ky,kz])


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
            kmult=sympy.matrices.Matrix.multiply(k,diffr.T)[0]
            t1=(ckdi*cmkdj*exp(-I*kmult)+cmkdi*ckdj*exp(I*kmult))/2
            t2=(cki*cmkj*exp(I*kmult)+cmki*ckj*exp(-I*kmult))/2
            t3=(ckdi*ckj*exp(-I*kmult)+cmkdi*cmkj*exp(I*kmult))/2
            t4=(cki*ckdj*exp(I*kmult)+cmki*cmkdj*exp(-I*kmult))/2
            t5=(ckdj*ckj+cmkdj*cmkj)/2
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
    return Hlin#.expand()        
    

    

if __name__=="__main__":
    N_atoms_uc=1
    Hdef=generate_hdef(N_atoms_uc)
    print 'Hdef',Hdef
    print Hdef.atoms(sympy.Symbol)
    B=sympy.Symbol('B')
    cd1=sympy.Symbol('cd1',commutative=False,real=True)
    c1=sympy.Symbol('c1',commutative=False,real=True)
    i=1
    j=1
    Hdef=cd1*c1*B
    cdilabel="cd%d"%(i,)
    cjlabel='c%d'%(j,)
    cdi=sympy.Symbol(cdilabel,commutative=False,real=True)
    cj=sympy.Symbol(cjlabel,commutative=False,real=True)
    print Hdef.subs(cd1*c1,B)
    print Hdef.subs(cdi*cj,B)

    