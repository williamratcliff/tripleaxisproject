import numpy as N

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
        D=sympy.Symbol('D',real=True)
        spin0=N.matrix([[1,0,0],[0,1,0],[0,0,1]],'float64')
        pos0=[0,0,0]
        neighbors=[1]
        interactions=[0]
        cell=0
        int_cell=[5,21]
        atom0=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=0,cell=cell,int_cell=int_cell,Dz=D)
        
        pos0=[1,0,0]
        spin0=N.matrix([[-1,0,0],[0,1,0],[0,0,-1]],'float64')
        neighbors=[0]
        interactions=[0]
        cell=5
        int_cell=[0]
        atom1=atom(spin=spin0,pos=pos0,neighbors=neighbors,interactions=interactions,label=1,cell=cell,int_cell=int_cell,Dz=D)
        
        atomlist=[atom0,atom1]


 
       
    return atomlist

def get_tokenized_line(myfile,returnline=['']):
        lineStr=myfile.readline()
        returnline[0]=lineStr.rstrip()
        strippedLine=lineStr.lower().rstrip()
        tokenized=strippedLine.split()

        return tokenized


def read_interactions(myfilestr):
    myfile = open(myfilestr, 'r')
    myFlag=True
        #self.metadata={}
    returnline=['']
    jmats=[]
    jnums=[]
    while myFlag:
        tokenized=get_tokenized_line(myfile,returnline=returnline)
        print tokenized
        if not(tokenized):
            break
        #print tokenized
        if tokenized==[]:
            break
        if tokenized[0]=='#number':
            while 1:
                tokenized=get_tokenized_line(myfile,returnline=returnline)
                print 'intoken ',tokenized
                if tokenized==[]:
                    break
                if tokenized[0]!='#atomnumber':
                    #print tokenized[0]
                    jnum=float(tokenized[0])
                    j11=float(tokenized[1])
                    j12=float(tokenized[2])
                    j13=float(tokenized[3])
                    j21=float(tokenized[4])
                    j22=float(tokenized[5])
                    j23=float(tokenized[6])
                    j31=float(tokenized[7])
                    j32=float(tokenized[8])
                    j33=float(tokenized[9]) 
                    jij=N.matrix([[j11,j12,j13],[j21,j22,j23],[j31,j32,j33]],'Float64')
                    jnums.append(jnum)
                    jmats.append(jij)
                else:
                    while 1:
                        tokenized=get_tokenized_line(myfile,returnline=returnline)
                        if not(tokenized):
                            break
                        atom_num=tokenized[0]
                        x,y,z=float(tokenized[1]),float(tokenized[2]),float(tokenized[3])
                        Dx,Dy,Dz=float(tokenized[4]),float(tokenized[5]),float(tokenized[6])
                        for i in range(7,1,len(tokenized)):
                            interacting_spin=float(tokenized[i])
                            interaction_matrix=float(tokenized[i+1])
                            
                                       
    myfile.close()
    
    
    
    
    
if __name__=="__main__":
    myfilestr=r'c:\montecarlo.txt'
    read_interactions(myfilestr)