import sympy
import numpy as N
I=sympy.I


class atom:
    def __init__(self,spin=[0,0,1],pos=[0,0,0],neighbors=[],interactions=[],label=0,Dx=0,Dy=0,Dz=0):
        self.spin=spin
        self.pos=N.array(pos)
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

def holstein(Hdef):
        S = sympy.Symbol('S')
        Hdef=Hdef.expand()
        Hdef=Hdef.as_poly(S)
        p = sympy.Wild('p',exclude='S')
        q = sympy.Wild('q',exclude='S')
        r = sympy.Wild('r',exclude='S')
        l = sympy.Wild('l',exclude='S')
        #Hlin=Hdef.match(p*S**2+S*q+l)
        #print Hdef
        #print Hdef.coeffs[0]
        #print Hdef.coeffs[1]
        Hlin=Hdef.coeffs[0]*S**2+Hdef.coeffs[1]*S
        return Hlin


def fouriertransform(atom_list,Jij,Hlin,k):
    N_atoms=len(atom_list)
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    for i in range(N_atoms):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False)
        ri=atom_list[i].pos
        for j in range(N_atoms):
            rj=atom_list[j].pos
            #rj=atom_list[atom_list[i].neighbors[j]].pos
            cj=sympy.Symbol("c%d"%(j,),commutative=False)
            cdj=sympy.Symbol("cd%d"%(j,),commutative=False)
            ckj=sympy.Symbol("ck%d"%(j,),commutative=False)
            ckdj=sympy.Symbol("ckd%d"%(j,),commutative=False)
            cmkj=sympy.Symbol("cmk%d"%(j,),commutative=False)
            cmkdj=sympy.Symbol("cmkd%d"%(j,),commutative=False)
            diffr=ri-rj
            kmult=N.dot(k,diffr)
            t1=1.0/2*(ckdi*cmkdj*sympy.exp(-I*kmult)+
                      cmkdi*ckdj*sympy.exp(I*kmult)               )
            t2=1.0/2*(cki*cmkj*sympy.exp(I*kmult)+
                      cmki*ckj*sympy.exp(-I*kmult)
                      )
            t3=1./2*(ckdi*ckj*sympy.exp(-I*kmult)+
                     cmkdi*cmkj*sympy.exp(I*kmult)
                     )
            t4=1./2*(cki*ckdj*sympy.exp(I*kmult)+
                     cmki*cmkdj*sympy.exp(-I*kmult)
                     )
            t5=1.0/2*(ckdj*ckj+cmkdj*cmkj
                      )
            print 'i',i,'j',j
            Hlin=Hlin.subs(cdi*cdj,t1)
            Hlin=Hlin.subs(ci*cj,t2)
            Hlin=Hlin.subs(cdi*cj,t3)
            Hlin=Hlin.subs(ci*cdj,t4)
            Hlin=Hlin.subs(cdj*cj,t5)
    #print t1
    return Hlin


def applycommutation(atom_list,Jij,Hfou,k):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""
    N_atoms=len(atom_list)
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    for i in range(N_atoms):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False)
        for j in range(N_atoms):
            cj=sympy.Symbol("c%d"%(j,),commutative=False)
            cdj=sympy.Symbol("cd%d"%(j,),commutative=False)
            ckj=sympy.Symbol("ck%d"%(j,),commutative=False)
            ckdj=sympy.Symbol("ckd%d"%(j,),commutative=False)
            cmkj=sympy.Symbol("cmk%d"%(j,),commutative=False)
            cmkdj=sympy.Symbol("cmkd%d"%(j,),commutative=False)
            if i==j:
                Hfou.subs(cki*ckdj,ckdj*cki+1)
                Hfou.subs(cki*ckdj,ckdj*cki)
                Hfou.subs(cmkdi*cmkj,cmkj*cmkdi+1)
                Hfou.subs(cmkdi*cmkj,cmkj*cmkdi)
            else:
                Hfou.subs(cki*cmkj,cmkj*cki)
                Hfou.subs(cmkdi*ckdj,ckdj*cmkdi)
    
    
    return Hfou

def gen_operator_table(atom_list):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""
    N_atoms=len(atom_list)
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    operator_table=[]
    operator_table_minus=[]
    
    for i in range(N_atoms):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False)
        operator_table.append(cki)
        operator_table_minus.append(cmkdi)
    operator_table=[operator_table,operator_table_minus]
    operator_table=N.ravel(operator_table)
    return operator_table


def gen_operator_table_dagger(atom_list):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""
    N_atoms=len(atom_list)
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    operator_table=[]
    operator_table_minus=[]
    
    for i in range(N_atoms):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False)
        operator_table.append(ckdi)
        operator_table_minus.append(cmki)
    operator_table=[operator_table,operator_table_minus]
    operator_table=N.ravel(operator_table)
    return operator_table

def gen_XdX(atom_list,operator_table,operator_table_dagger,Hcomm):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""
    N_atoms=len(atom_list)
    #Hdef=0
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    exclude_list=[]
    coeff_list=[]
    for i in range(2*N_atoms):
        for j in range(2*N_atoms):
            coeff=operator_table_dagger[i]*operator_table[j]
            #print 'coeff',coeff,'i',i,'j',j
            exclude_list.append(coeff)
    for i in range(len(exclude_list)):
        coeff_list.append(sympy.Wild('A%d'%(i),exclude=exclude_list))
    coeff_list=N.array(coeff_list)
    exclude_list2=N.array(exclude_list)
    expr=N.sum(coeff_list*exclude_list2)
    print expr
    #B=sympy.Wild('B',exclude=exclude_list)
    Hcomm=Hcomm.expand()
    print Hcomm
    resultdict=Hcomm.match(expr)
    print 'done'
    print 'res',resultdict
    XdX=sympy.Matrix((3,3))       
    return XdX


    
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
        N_atoms=2
        Sabn=generate_sabn(N_atoms)
        print Sabn[0]
        Sxyz=generate_sxyz(N_atoms)
        atom_list=generate_atoms()
        Jij=[N.matrix([[1,0,0],[0,1,0],[0,0,1]])]
        Hdef=generate_hdef(atom_list,Jij,Sabn)
        print 'Hdef',Hdef
        Hlin=holstein(Hdef)
        print 'Hlin',Hlin
        kx=sympy.Symbol('kx')
        ky=sympy.Symbol('ky')
        kz=sympy.Symbol('kz')
        k=[kx,ky,kz]
        Hfou=fouriertransform(atom_list,Jij,Hlin,k)
        print 'Hfou',Hfou
        Hcomm=applycommutation(atom_list,Jij,Hfou,k)
        print 'Hcomm',Hcomm
        operator_table=gen_operator_table(atom_list)
        print 'optable',operator_table
        operator_table_dagger=gen_operator_table_dagger(atom_list)
        print 'optable_dagger',operator_table_dagger
        XdX=gen_XdX(atom_list,operator_table,operator_table_dagger,Hcomm)
        print 'XdX',XdX
        