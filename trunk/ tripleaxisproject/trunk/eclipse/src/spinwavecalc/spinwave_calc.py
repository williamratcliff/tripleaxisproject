import sympy
import numpy as N
#I=sympy.I
I=1.0j
pi=sympy.pi
from sympy import exp

#translations=[[0,0,0],
#              [0,0,1],[0,0,-1]
#              [0,1,0],[0,-1,0]
#              [0,1,1],[0,1,-1],[0,-1,1],[0,-1,-1]
#              [1,0,0],[-1,0,0],[]              
#              ]



def coeff(expr, term):
   expr = sympy.collect(expr, term)
   #print 'expr',expr
   symbols = list(term.atoms(sympy.Symbol))
   #print 'symbols',symbols
   w = sympy.Wild("coeff", exclude=symbols)
   #print 'w',w
   m = expr.match(w*term+sympy.Wild("rest"))
   #print 'm',m
   m2=expr.match(w*term)
   #print 'm2',m2
   res=False
   if m2!=None:
       #print 'm2[w]',m2[w]
       res=m2[w]*term==expr
   if m and res!=True:
       return m[w]
   #added the next two lines
   elif m2:
       return m2[w]
       
   
class atom:
    def __init__(self,spin=[0,0,1],pos=[0,0,0],neighbors=[],interactions=[],label=0,Dx=0,Dy=0,Dz=0,cell=0,int_cell=[]):
        self.spin=spin
        self.pos=N.array(pos)
        self.neighbors=neighbors
        self.interactions=interactions
        self.label=label
        self.Dx=Dx
        self.Dy=Dy
        self.Dz=Dz
        self.cell=cell
        self.int_cell=[]


def generate_atoms():



    if 1:
        spin0=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'float64')
        pos0=[0,0,0]
        neighbors=[1]
        interactions=[0]
        cell=0
        int_cell=[5,21]
        atom0=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=0,cell=cell,int_cell=int_cell)
        
        pos0=[1,0,0]
        spin0=N.matrix([[1,0,0],[0,-1,0],[0,0,-1]],'float64')
        neighbors=[0]
        interactions=[0]
        cell=5
        int_cell=[0]
        atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1,cell=cell,int_cell=int_cell)
        
        atomlist=[atom0,atom1]


 
       
    return atomlist



def generate_sabn(N_atoms):
    Sabn=[]
    S=sympy.Symbol("S",real=True)
    for i in range(N_atoms):
        c=sympy.Symbol("c%d"%(i,),commutative=False)
        cd=sympy.Symbol("cd%d"%(i,),commutative=False)
        curr=[sympy.sqrt(S/2)*(c+cd),sympy.sqrt(S/2)*(c-cd)/I,S-cd*c]
        Sabn.append(curr)
    return Sabn


def generate_sxyz(Sabn,atomlist):
    Sxyz=[]
    i=0
    for currS in Sabn:
        print 'Currs', currS
        print 'currspin', atomlist[i].spin
        currS_transpose=N.reshape(currS,(3,1))
        tempS=atomlist[i].spin*currS_transpose
        tempS=N.array(tempS)
        tempS=N.ravel(tempS)
        Sxyz.append(tempS)
        i=i+1
    return Sxyz




def generate_translations():
    translations=[[0,0,0]]
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                if i==0 and j==0 and k==0:
                    continue
                else:
                    translations.append([i,j,k])
    return translations

def generate_sabnt(N_atoms,t=''):
    Sabn=[]
    S=sympy.Symbol("S",real=True)
    for i in range(N_atoms):
        c=sympy.Symbol("c%d%s"%(i,t),commutative=False)
        cd=sympy.Symbol("cd%d%s"%(i,t),commutative=False)
        curr=[sympy.sqrt(S/2)*(c+cd),sympy.sqrt(S/2)*(c-cd)/I,S-cd*c]
        Sabn.append(curr)
    return Sabn





def generate_hdef(atom_list,Jij,Sxyz,translations,N_atoms_uc,N_atoms):
    N_atoms=len(atom_list)
    Hdef=0
    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        for j in range(N_int):
            Hij=Sxyz[i]*Jij[atom_list[i].interactions[j]]#
            Sxyz_transpose=N.matrix(Sxyz[atom_list[i].neighbors[j]])
            Sxyz_transpose=N.reshape(Sxyz_transpose,(3,1))
            Hij=Hij*Sxyz_transpose
            Hij=-Hij-atom_list[i].Dx*Sxyz[i][0]**2-atom_list[i].Dy*Sxyz[i][1]**2-atom_list[i].Dz*Sxyz[i][2]**2
            Hdef=Hdef+Hij
    return Hdef[0][0]


def holstein(Hdef):
        S = sympy.Symbol('S',real=True)
        Hdef=Hdef.expand()
        #Hdef=Hdef.as_poly(S)
        p = sympy.Wild('p',exclude='S')
        q = sympy.Wild('q',exclude='S')
        r = sympy.Wild('r',exclude='S')
        l = sympy.Wild('l',exclude='S')
        #Hlin=Hdef.coeffs[0]*S**2+Hdef.coeffs[1]*S
        S2coeff=coeff(Hdef,S**2)
        Scoeff=coeff(Hdef,S)
        Hlin=None
        #Hlin=coeff(Hdef,S**2)*S**2+coeff(Hdef,S)*S
        #print 'S2Coeff', S2coeff
        #print 'Scoeff',Scoeff
        if Scoeff!=None and S2coeff!=None:
            Hlin=coeff(Hdef,S**2)*S**2+coeff(Hdef,S)*S
        elif Scoeff==None and S2coeff!=None:
            Hlin=coeff(Hdef,S**2)*S**2
        elif Scoeff!=None and S2coeff==None:
            #print 'S'
            Hlin=coeff(Hdef,S)*S
        return Hlin


def fouriertransform(atom_list,Jij,Hlin,k,N_atoms_uc,N_atoms):
    #N_atoms=len(atom_list)
    #N_atoms_uc=1
    #N_atoms_uc=N_atoms
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False,real=True)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False,real=True)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False,real=True)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False,real=True)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False,real=True)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False,real=True)
        ri=atom_list[i].pos
        N_int=len(atom_list[i].interactions)
        for j in range(N_atoms):
            rj=atom_list[j].pos
            j2=i#atom_list[i].neighbors[j]
            cj=sympy.Symbol("c%d"%(j,),commutative=False,real=True)
            cdj=sympy.Symbol("cd%d"%(j,),commutative=False,real=True)
            ckj=sympy.Symbol("ck%d"%(j2,),commutative=False,real=True)
            ckdj=sympy.Symbol("ckd%d"%(j2,),commutative=False,real=True)
            cmkj=sympy.Symbol("cmk%d"%(j2,),commutative=False,real=True)
            cmkdj=sympy.Symbol("cmkd%d"%(j2,),commutative=False,real=True)
            diffr=ri-rj
            kmult=N.dot(k,diffr)
            t1=1.0/2*(ckdi*cmkdj*exp(-1j*kmult)+
                      cmkdi*ckdj*exp(1j*kmult)               )
            t2=1.0/2*(cki*cmkj*exp(1j*kmult)+
                      cmki*ckj*exp(-1j*kmult)
                      )
            t3=1./2*(ckdi*ckj*exp(-1j*kmult)+
                     cmkdi*cmkj*exp(1j*kmult)
                     )
            t4=1./2*(cki*ckdj*exp(1j*kmult)+
                     cmki*cmkdj*exp(-1j*kmult)
                     )
            t5=1.0/2*(ckdj*ckj+cmkdj*cmkj
                      )
            #print 'i',i,'j',j
            Hlin=Hlin.subs(cdi*cdj,t1)
            Hlin=Hlin.subs(ci*cj,t2)
            Hlin=Hlin.subs(cdi*cj,t3)
            Hlin=Hlin.subs(ci*cdj,t4)
            Hlin=Hlin.subs(cdj*cj,t5)
    #print t1
    return Hlin.expand()


def applycommutation(atom_list,Jij,Hfou,k,N_atoms_uc,N_atoms):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""

    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False,real=True)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False,real=True)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False,real=True)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False,real=True)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False,real=True)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False,real=True)
        for j in range(N_atoms_uc):
            cj=sympy.Symbol("c%d"%(j,),commutative=False,real=True)
            cdj=sympy.Symbol("cd%d"%(j,),commutative=False,real=True)
            ckj=sympy.Symbol("ck%d"%(j,),commutative=False,real=True)
            ckdj=sympy.Symbol("ckd%d"%(j,),commutative=False,real=True)
            cmkj=sympy.Symbol("cmk%d"%(j,),commutative=False,real=True)
            cmkdj=sympy.Symbol("cmkd%d"%(j,),commutative=False,real=True)
            if i==j:
                Hfou=Hfou.subs(cki*ckdj,ckdj*cki+1)
                
                Hfou=Hfou.subs(cmkdi*cmkj,cmkj*cmkdi+1)
                
            else:
                Hfou=Hfou.subs(cki*ckdj,ckdj*cki)
                Hfou=Hfou.subs(cmkdi*cmkj,cmkj*cmkdi)
            Hfou=Hfou.subs(cki*cmkj,cmkj*cki)
                
            Hfou=Hfou.subs(cmkdi*ckdj,ckdj*cmkdi)
    
    
    return Hfou.expand()

def gen_operator_table(atom_list,N_atoms_uc):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""

    operator_table=[]
    operator_table_minus=[]
    
    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False,real=True)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False,real=True)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False,real=True)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False,real=True)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False,real=True)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False,real=True)
        operator_table.append(cki)
        operator_table_minus.append(cmkdi)
    operator_table=[operator_table,operator_table_minus]
    operator_table=N.ravel(operator_table)
    return operator_table


def gen_operator_table_dagger(atom_list,N_atoms_uc):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""

    operator_table=[]
    operator_table_minus=[]
    
    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False,real=True)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False,real=True)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False,real=True)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False,real=True)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False,real=True)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False,real=True)
        operator_table.append(ckdi)
        operator_table_minus.append(cmki)
    operator_table=[operator_table,operator_table_minus]
    operator_table=N.ravel(operator_table)
    return operator_table

def gen_XdX(atom_list,operator_table,operator_table_dagger,Hcomm,N_atoms_uc):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""

    exclude_list=[]
    coeff_list=[]
    Hcomm=Hcomm.expand()
    XdX=sympy.zeros(2*N_atoms_uc)
    g=sympy.zeros(2*N_atoms_uc)
    for i in range(2*N_atoms_uc):
        curr_row=[]
        for j in range(2*N_atoms_uc):
            mycoeff=operator_table_dagger[i]*operator_table[j]
            exclude_list.append(mycoeff)
            currcoeff=coeff(Hcomm,mycoeff)
            if currcoeff!=None:
                XdX[i,j]=currcoeff
                curr_row.append(currcoeff)
            if i!=j:
                g[i,j]=0
            else:
                if i>=N_atoms_uc:
                    g[i,j]=-1
                else:
                    g[i,j]=1
 
    return XdX,g



def calculate_dispersion(atom_list,N_atoms_uc,N_atoms):
    Sabn=generate_sabn(N_atoms)        
    Sxyz=generate_sxyz(Sabn,atom_list)
    print 'Sabn',Sabn
    print 'Sxyz', Sxyz
        
    if 1:
        print len(translations)   
        J=sympy.Symbol('J',real=True)
        Jij=[N.matrix([[J,0,0],[0,J,0],[0,0,J]])]
        Hdef=generate_hdef(atom_list,Jij,Sxyz,translations,N_atoms_uc,N_atoms)
        print 'Hdef',Hdef
    #if 0:
        Hlin=holstein(Hdef)
        print 'Hlin',Hlin
        kx=sympy.Symbol('kx',real=True)
        ky=sympy.Symbol('ky',real=True)
        kz=sympy.Symbol('kz',real=True)
        k=[kx,ky,kz]
        Hfou=fouriertransform(atom_list,Jij,Hlin,k,N_atoms_uc,N_atoms)
        print 'Hfou',Hfou
    if 1:
        Hcomm=applycommutation(atom_list,Jij,Hfou,k,N_atoms_uc,N_atoms)
        print 'Hcomm',Hcomm
        operator_table=gen_operator_table(atom_list,N_atoms_uc)
        print 'optable',operator_table
        operator_table_dagger=gen_operator_table_dagger(atom_list,N_atoms_uc)
        print 'optable_dagger',operator_table_dagger
        XdX,g=gen_XdX(atom_list,operator_table,operator_table_dagger,Hcomm,N_atoms_uc)
        print 'XdX',XdX
        print 'g',g
        TwogH2=2*g*XdX
        print 'TwogH2',TwogH2
        if 1:
            #print 'calculating'
            x=sympy.Symbol('x')
            #eigspoly=TwogH2.berkowitz_charpoly(x)
            #print 'eigspoly'
            #print 'eigs poly',eigspoly
            print 'recalculating'
            eigs=TwogH2.eigenvals()
            #x=sympy.Symbol('x')
            #eigs=TwogH2.berkowitz_charpoly(x)
            print 'eigs', eigs
            keys=eigs.keys()
            print 'key',keys[0]
            #print keys[0].expand(complex=True)
            #print TwogH2.charpoly(x)
            eigs=TwogH2.eigenvals()
            #print 'eigenvalues', sympy.simplify(eigs[1][0])        
        S=sympy.Symbol('S',real=True)
        TwogH2=TwogH2.subs(J,1.0)
        TwogH2=TwogH2.subs(S,1.0)
        TwogH2=TwogH2.subs(kx,2*pi/4)
        #I=sympy.Symbol('I')
        TwogH2=TwogH2.subs(I,1.0j)
        Ntwo=N.array(TwogH2)
        print 'Ntwo'
        print Ntwo
        if 1:
            import scipy.linalg
            l,v=scipy.linalg.eig(Ntwo)
            print l[0]
        #print N.linalg.eigvals(Ntwo)
        #eigs=TwogH2.eigenvals()
        #print 'eigs', eigs
        #print 'eigenvalues', sympy.simplify(eigs[1][0])
        return eigs 



    
if __name__=='__main__':
    if 1:
        translations=generate_translations()
        atom_list=generate_atoms()
        N_atoms_uc=1
        N_atoms=2
        calculate_dispersion(atom_list,N_atoms_uc,N_atoms)
