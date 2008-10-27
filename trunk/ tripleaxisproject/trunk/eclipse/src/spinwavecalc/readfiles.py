

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
    while myFlag:
        tokenized=get_tokenized_line(myfile,returnline=returnline)
        print tokenized
        if not(tokenized):
            break
        #print tokenized
        if tokenized==[]:
            tokenized=['']
        if tokenized[0].lower()=="#Date".lower():
            pass
    myfile.close()
    
    
    
    
    
if __name__=="__main__":
    myfilestr=r'c:\montecarlo.txt'
    read_interactions(myfilestr)