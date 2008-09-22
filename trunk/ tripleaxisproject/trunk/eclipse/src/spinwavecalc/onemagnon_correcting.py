import sympy
import numpy as N
I=sympy.I
pi=N.pi

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

def generate_sabn(N_atoms):
    Sabn=[]
    S=sympy.Symbol("S",real=True)
    for i in range(N_atoms):
        c=sympy.Symbol("c%d"%(i,),commutative=False)
        cd=sympy.Symbol("cd%d"%(i,),commutative=False)
        curr=[sympy.sqrt(S/2)*(c+cd),sympy.sqrt(S/2)*(c-cd)/sympy.I,S-cd*c]
        Sabn.append(curr)
    return Sabn


def generate_sxyz(Sabn,atomlist):
    Sxyz=[]
    i=0
    for currS in Sabn:
        print 'Currs', currS
        print 'currspin', atomlist[i].spin
        currS_transpose=N.reshape(currS,(3,1))
        tempS=N.dot(atomlist[i].spin,currS_transpose)
        #tempS=N.reshape(tempS,(1,3))
        tempS=N.array(tempS)
        tempS=N.ravel(tempS)
        Sxyz.append(tempS)
        print 'tempS', tempS
    return Sxyz


#def generate_sxyz(N_atoms):
#    Sxyz=[]
#    for i in range(N_atoms):
#        S=sympy.Symbol("Sxyz%d"%(i,),commutative=False)
#        Sxyz.append(S)
#    return Sxyz





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
        curr=[sympy.sqrt(S/2)*(c+cd),sympy.sqrt(S/2)*(c-cd)/sympy.I,S-cd*c]
        Sabn.append(curr)
    return Sabn





def generate_hdef(atom_list,Jij,Sxyz,translations,N_atoms_uc,N_atoms):
    N_atoms=len(atom_list)
    Hdef=0
    #print 'atom_list',atom_list
    print 'Sxyz',Sxyz
    #print 'Jij',Jij
    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        for j in range(N_int):
            #print Jij[atom_list[i].interactions[j]]
            #print '1',Sxyz[i]
            #print Sxyz[atom_list[i].neighbors[j]]
            #print 'i',i,'j',j
            #currj=Jij[atom_list[i].interactions[j]]
            #print Jij[0]
            Hij=Sxyz[i]*Jij[atom_list[i].interactions[j]]#
            #print type(Hij)
            Sxyz_transpose=N.matrix(Sxyz[atom_list[i].neighbors[j]])
            Sxyz_transpose=N.reshape(Sxyz_transpose,(3,1))
            Hij=Hij*Sxyz_transpose
            Hij=-Hij-atom_list[i].Dx*Sxyz[i][0]**2-atom_list[i].Dy*Sxyz[i][1]**2-atom_list[i].Dz*Sxyz[i][2]**2
            Hdef=Hdef+Hij
    return Hdef[0][0]

def generate_atoms():
    if 0:
        spin0=N.matrix([[0,0,1],[0,1,0],[0,0,1]],'float64')
        pos0=[0,0,0]
        neighbors=[1,2]
        interactions=[0,0]
        cell=0
        int_cell=[5,21]
        atom0=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=0,cell=cell,int_cell=int_cell)
        
        pos0=[-1,0,0]
        neighbors=[0]
        interactions=[0]
        cell=5
        int_cell=[0]
        atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1,cell=cell,int_cell=int_cell)
        
        pos0=[1,0,0]
        neighbors=[0]
        interactions=[0]
        cell=22
        int_cell=[0]
        atom2=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=2,cell=cell,int_cell=int_cell)
        atomlist=[atom0,atom1,atom2]
        
    if 1:
        spin0=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'float64')
        pos0=[0,0,0]
        neighbors=[1,2]
        interactions=[0,0]
        cell=0
        int_cell=[5,21]
        atom0=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=0,cell=cell,int_cell=int_cell)
        
        pos0=[0.5,0,0]
        spin0=N.matrix([[-1,0,0],[0,-1,0],[0,0,-1]],'float64')
        neighbors=[0,3]
        interactions=[0,0]
        cell=5
        int_cell=[0]
        atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1,cell=cell,int_cell=int_cell)
        
        pos0=[-.5,0,0]
        spin0=N.matrix([[-1,0,0],[0,-1,0],[0,0,-1]],'float64')
        neighbors=[0]
        interactions=[0]
        cell=5
        int_cell=[0]
        atom2=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1,cell=cell,int_cell=int_cell)

        pos0=[1,0,0]
        spin0=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'float64')
        neighbors=[1]
        interactions=[0]
        cell=5
        int_cell=[0]
        atom3=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1,cell=cell,int_cell=int_cell)

        atomlist=[atom0,atom1,atom2,atom3]
        
    if 0:
        spin0=[0,0,1]
        pos0=[0,0,0]
        neighbors=[1]
        interactions=[0]
        atom0=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=0)
        
        pos0=[1,0,0]
        neighbors=[0]
        interactions=[0]
        atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1)
        atomlist=[atom0,atom1]
       
    return atomlist

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
        ci=sympy.Symbol("c%d"%(i,),commutative=False)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False)
        ri=atom_list[i].pos
        N_int=len(atom_list[i].interactions)
        for j in range(N_atoms):
            rj=atom_list[j].pos
            j2=i#atom_list[i].neighbors[j]
            #rj=atom_list[atom_list[i].neighbors[j]].pos
            #print 'j2',j2
            #rj=atom_list[atom_list[i].neighbors[j]].pos
            cj=sympy.Symbol("c%d"%(j,),commutative=False)
            cdj=sympy.Symbol("cd%d"%(j,),commutative=False)
            ckj=sympy.Symbol("ck%d"%(j2,),commutative=False)
            ckdj=sympy.Symbol("ckd%d"%(j2,),commutative=False)
            cmkj=sympy.Symbol("cmk%d"%(j2,),commutative=False)
            cmkdj=sympy.Symbol("cmkd%d"%(j2,),commutative=False)
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
    #N_atoms=len(atom_list)
    #N_atoms_uc=1
    #N_atoms_uc=N_atoms
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    for i in range(N_atoms_uc):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(i,),commutative=False)
        cdi=sympy.Symbol("cd%d"%(i,),commutative=False)
        cki=sympy.Symbol("ck%d"%(i,),commutative=False)
        ckdi=sympy.Symbol("ckd%d"%(i,),commutative=False)
        cmki=sympy.Symbol("cmk%d"%(i,),commutative=False)
        cmkdi=sympy.Symbol("cmkd%d"%(i,),commutative=False)
        for j in range(N_atoms_uc):
            cj=sympy.Symbol("c%d"%(j,),commutative=False)
            cdj=sympy.Symbol("cd%d"%(j,),commutative=False)
            ckj=sympy.Symbol("ck%d"%(j,),commutative=False)
            ckdj=sympy.Symbol("ckd%d"%(j,),commutative=False)
            cmkj=sympy.Symbol("cmk%d"%(j,),commutative=False)
            cmkdj=sympy.Symbol("cmkd%d"%(j,),commutative=False)
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
    #N_atoms=len(atom_list)
    #N_atoms_uc=1
    #N_atoms_uc=N_atoms
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    operator_table=[]
    operator_table_minus=[]
    
    for i in range(N_atoms_uc):
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


def gen_operator_table_dagger(atom_list,N_atoms_uc):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""
    #N_atoms=len(atom_list)
    #N_atoms_uc=1
    #N_atoms_uc=N_atoms
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    operator_table=[]
    operator_table_minus=[]
    
    for i in range(N_atoms_uc):
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

def gen_XdX(atom_list,operator_table,operator_table_dagger,Hcomm,N_atoms_uc):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""
    #N_atoms=len(atom_list)
    #N_atoms_uc=1
    #N_atoms_uc=N_atoms
    #Hdef=0
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    exclude_list=[]
    coeff_list=[]
    Hcomm=Hcomm.expand()
    XdX=sympy.zeros(2*N_atoms_uc)
    g=sympy.zeros(2*N_atoms_uc)
    #print 'XdX',XdX    
    for i in range(2*N_atoms_uc):
        curr_row=[]
        for j in range(2*N_atoms_uc):
            mycoeff=operator_table_dagger[i]*operator_table[j]
            #print 'coeff',coeff,'i',i,'j',j
            exclude_list.append(mycoeff)
            currcoeff=coeff(Hcomm,mycoeff)
            #print 'found',currcoeff
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
                
            #else:
            #    XdX[i,j]=0
            #    curr_row.append(0)
        #XdX[i]=curr_row
                
            
    #for i in range(len(exclude_list)):
    #    coeff_list.append(sympy.Wild('A%d'%(i),exclude=exclude_list))
    #coeff_list2=N.array(coeff_list)
    #exclude_list2=N.array(exclude_list)
    #expr=N.sum(coeff_list2*exclude_list2)


    #resultdict=Hcomm.match(expr)
    #print 'done'
    #print 'res',resultdict
     
    #XdX[0][0]=1  
    return XdX,g



def multiply_ab(atom_list,Sxyz,a=0,b=0):
    N_atoms=len(atom_list)
    Sdef=0
    print 'atom_list',atom_list
    print 'Sxyz',Sxyz
    T=sympy.Symbol('T',commutative=False)
    Sij0=Sxyz[0][a]
    t=''
    c=sympy.Symbol("c%d%s"%(0,t),commutative=False)
    cd=sympy.Symbol("cd%d%s"%(0,t),commutative=False)
    t='t'
    ct=sympy.Symbol("c%d%s"%(0,t),commutative=False)
    cdt=sympy.Symbol("cd%d%s"%(0,t),commutative=False)
    Sij0=Sij0.subs(ct,c)
    Sij0=Sij0.subs(cdt,cd)  
    for i in range(1,N_atoms):
        Sdef=Sdef+Sij0*Sxyz[i][b]
 
    return Sdef




def Sfouriertransform(atom_list,Slin,k):
    N_atoms=len(atom_list)
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    kxp=sympy.Symbol('kxp')
    kyp=sympy.Symbol('kyp')
    kzp=sympy.Symbol('kzp')
    kp=[kxp,kyp,kzp]
    wk=sympy.Symbol("wk")
    wkp=sympy.Symbol("wkp")
    t=sympy.Symbol("t")
    ri=atom_list[0].pos
    
    kmult=N.dot(k,ri)
    #kmultp=N.dot(kp,ri)
    
    ci=sympy.Symbol("c%d"%(0,),commutative=False)
    cdi=sympy.Symbol("cd%d"%(0,),commutative=False)
    cki=sympy.Symbol("ck%d"%(0,),commutative=False)
    ckdi=sympy.Symbol("ckd%d"%(0,),commutative=False)
    
    t1=sympy.exp(I*kmult)*cki
    t3=sympy.exp(-I*(kmult))*ckdi
    Slin=Slin.subs(ci,t1)
    Slin=Slin.subs(cdi,t3)


    for i in range(1,N_atoms):
        N_int=len(atom_list[i].interactions)
        #ci=sympy.Symbol("c%d"%(i,),commutative=False)
        #cdi=sympy.Symbol("cd%d"%(i,),commutative=False)
        cit=sympy.Symbol("c%dt"%(i,),commutative=False)
        cdit=sympy.Symbol("cd%dt"%(i,),commutative=False)

        cki=sympy.Symbol("ck%d"%(0,),commutative=False)
        ckdi=sympy.Symbol("ckd%d"%(0,),commutative=False)
        
        ri=atom_list[i].pos
        kmult=N.dot(k,ri)
        kmultp=N.dot(kp,ri)
        
        t2=sympy.exp(I*(kmultp-wk*t))*cki
                     
        
        t4=sympy.exp(-I*(kmultp-wkp*t))*ckdi
        
        Slin=Slin.subs(cit,t2)
        
        Slin=Slin.subs(cdit,t4)
    
    #Note that I have already assumed that k=kp because I didn't include cqp terms
    Slin=Slin.expand()
    Slin=Slin.subs(wkp,wk)
    Slin=Slin.subs(kxp,kx)
    Slin=Slin.subs(kyp,ky)
    Slin=Slin.subs(kzp,kz)
    return Slin

def Sapplycommutation(atom_list,Sfou,k):
    """Operate commutation relations to put all the 2nd order term as ckd**ck, cmk**cmkd, cmk**ck and ckd**cmkd form"""
    N_atoms=len(atom_list)
    #Hdef=0
    #print 'atom_list',atom_list
    #print 'Sxyz',Sxyz
    #print 'Jij',Jij
    for i in range(N_atoms):
        N_int=len(atom_list[i].interactions)
        ci=sympy.Symbol("c%d"%(0,),commutative=False)
        cdi=sympy.Symbol("cd%d"%(0,),commutative=False)
        cki=sympy.Symbol("ck%d"%(0,),commutative=False)
        ckdi=sympy.Symbol("ckd%d"%(0,),commutative=False)
        cmki=sympy.Symbol("cmk%d"%(0,),commutative=False)
        cmkdi=sympy.Symbol("cmkd%d"%(0,),commutative=False)
        for j in range(N_atoms):
            cj=sympy.Symbol("c%d"%(0,),commutative=False)
            cdj=sympy.Symbol("cd%d"%(0,),commutative=False)
            ckj=sympy.Symbol("ck%d"%(0,),commutative=False)
            ckdj=sympy.Symbol("ckd%d"%(0,),commutative=False)
            cmkj=sympy.Symbol("cmk%d"%(0,),commutative=False)
            cmkdj=sympy.Symbol("cmkd%d"%(0,),commutative=False)
            Sfou=Sfou.expand()
            if i==j:
                Sfou=Sfou.subs(cki*ckdj,ckdj*cki+1)
                #Sfou.subs(cmkdi*cmkj,cmkj*cmkdi)
                nkj=sympy.Symbol("nk%d"%(j,),commutative=True)
                
            else:
                Sfou=Sfou.subs(cki*ckdj,ckdj*cki) #just added
            
            Sfou=Sfou.expand()
            nkj=sympy.Symbol("nk%d"%(j,),commutative=True)
            Sfou=Sfou.subs(ckdj*cki,nkj)
            Sfou=Sfou.expand()
            Sfou=Sfou.subs(cki*ckj,0)
            Sfou=Sfou.expand()
            Sfou=Sfou.subs(ckdi*ckdj,0)

    
    
    return Sfou
    
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
        translations=generate_translations()
        print 'translations',translations
        print translations[0]
        print translations[5]
        print translations[22]
        atom_list=generate_atoms()
        N_atoms_uc=2
        N_atoms=4
        Sabn=generate_sabn(N_atoms)
        print Sabn[0]
        Sxyz=generate_sxyz(Sabn,atom_list)
        print len(translations)   
        J=sympy.Symbol('J',real=True)
        Jij=[N.matrix([[J,0,0],[0,J,0],[0,0,J]])]
        #Hdef=generate_hdef(atom_list,Jij,Sabn,translations,N_atoms_uc,N_atoms)
        Hdef=generate_hdef(atom_list,Jij,Sxyz,translations,N_atoms_uc,N_atoms)
        print 'Hdef',Hdef
    #if 0:
        Hlin=holstein(Hdef)
        print 'Hlin',Hlin
        kx=sympy.Symbol('kx')
        ky=sympy.Symbol('ky')
        kz=sympy.Symbol('kz')
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
            print 'calculating'
            x=sympy.Symbol('x')
            eigspoly=TwogH2.berkowitz_charpoly(x)
            print 'eigspoly'
            print 'eigs poly',eigspoly
            print 'recalculating'
            eigs=TwogH2.eigenvals()
            #x=sympy.Symbol('x')
            #eigs=TwogH2.berkowitz_charpoly(x)
            print 'eigs', eigs
            #print TwogH2.charpoly(x)
            eigs=TwogH2.eigenvals()
            #print 'eigenvalues', sympy.simplify(eigs[1][0])        
        S=sympy.Symbol('S',real=True)
        TwogH2=TwogH2.subs(J,1.0)
        TwogH2=TwogH2.subs(S,1.0)
        TwogH2=TwogH2.subs(kx,2*pi)
        Ntwo=N.array(TwogH2)
        #print Ntwo
        if 0:
            import scipy.linalg
            l,v=scipy.linalg.eig(Ntwo)
            print l  
        #print N.linalg.eigvals(Ntwo)
        #eigs=TwogH2.eigenvals()
        #print 'eigs', eigs
        #print 'eigenvalues', sympy.simplify(eigs[1][0])
    if 0:
        print 'one magnon'
        print ''
        print ''
        Sabnt=generate_sabnt(N_atoms,t='t')
        SzSz=sympy.expand(multiply_ab(atom_list,Sabnt,a=2,b=2))
        print 'mult SzSz',SzSz
        SxSx=sympy.expand(multiply_ab(atom_list,Sabnt,a=0,b=0))
        SxSy=sympy.expand(multiply_ab(atom_list,Sabnt,a=0,b=1))
        print 'mult SxSx',SxSx
        SzSz_lin=holstein(sympy.expand(SzSz))
        print 'lin zz',SzSz_lin
        SxSx_lin=holstein(sympy.expand(SxSx))
        SxSy_lin=holstein(sympy.expand(SxSy))
        print 'lin xy',SzSz_lin
        print 'lin xx',SxSx_lin        
        if SzSz_lin!=None:
            SzSz_fou=Sfouriertransform(atom_list,SzSz_lin,k)
            print 'fourier', SzSz_fou
            Scomm=Sapplycommutation(atom_list,SzSz_fou,k)
            print 'Scomm',Scomm
            SxSx_fou=Sfouriertransform(atom_list,SxSx_lin,k)
            print 'fourier x',SxSx_fou
            Scommx=Sapplycommutation(atom_list,SxSx_fou,k)
            print 'Scommx',sympy.simplify(Scommx)
            SxSy_fou=Sfouriertransform(atom_list,SxSy_lin,k)
            print 'fourier xy',SxSy_fou
            Scommxy=Sapplycommutation(atom_list,SxSy_fou,k)
            print 'Scommxy',Scommxy            