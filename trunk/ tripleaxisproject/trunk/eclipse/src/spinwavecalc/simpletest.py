import sympy

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


if __name__=="__main__":
    W=sympy.Symbol('W',real=True,commutative=True)
    X=sympy.Symbol('X',real=True,commutative=True)
    P=sympy.Symbol('P',real=True,commutative=True)
    a=sympy.Symbol('a',real=True,commutative=False)
    b=sympy.Symbol('b',real=True,commutative=False)
    c=sympy.Symbol('c',real=True,commutative=False)
    d=sympy.Symbol('d',real=True,commutative=False)
    e=sympy.Symbol('e',real=True,commutative=False)
    f=sympy.Symbol('f',real=True,commutative=False)
    x=sympy.Symbol('x',real=True,commutative=False)
    spin0=sympy.matrices.Matrix([[W,-X,0],[W,X,0],[0,0,1]])
    spin2=sympy.matrices.Matrix([[1,0,0],[0,1,0],[0,0,1]])
    Jij=[sympy.matrices.Matrix([[1,0,0],[0,1,0],[0,0,1]])]
    #Hij=Sxyz[i]*Jij[atom_list[i].interactions[j]]
    Sabn1=sympy.matrices.Matrix([a,b,c]).reshape(3,1)
    Sxyz1=spin0*Sabn1
    Sabn2=sympy.matrices.Matrix([d,e,f]).reshape(3,1)
    Sxyz2=spin2*Sabn2
    #SxyzT=Sxyz.reshape(1,3)
    Hij=Sxyz1.T*Jij[0]*P
    print Hij
    Hij=Hij*Sxyz2
    print Hij
    Hij=Hij[0]
    print Hij.subs(W*a,x)