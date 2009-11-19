import numpy as np
from enthought.mayavi import mlab

def gen_as():
    xas=0.361
    pos_at=np.array([[0.0000,   0.0000,   1-xas],\
                [0.0000,   0.0000,   xas],\
                [0.0000,   0.5000,   .5-xas],
                [0.0000,   0.5000,   .5+xas],\
                [0.5000,   0.0000,   .5-xas],\
                [0.5000,   0.0000,   .5+xas],\
                [0.5000,   0.5000,   1-xas],\
                [0.5000,   0.5000,   xas]])      
    x=pos_at[:,0]
    y=pos_at[:,1]
    z=pos_at[:,2]
    return x,y,z

def gen_fe():
    pos_at=np.array([[0.2500,   0.2500,   0.2500],\
                [0.7500,   0.7500,   0.2500],\
                [0.7500,   0.2500,   0.7500],
                [0.2500,   0.7500,   0.7500],\
                [0.7500,   0.7500,   0.7500],\
                [0.2500,   0.2500,   0.7500],\
                [0.2500,   0.7500,   0.2500],\
                [0.7500,   0.2500,   0.2500]])      
    x=pos_at[:,0]
    y=pos_at[:,1]
    z=pos_at[:,2]
    return x,y,z


def gen_sr():
    pos_at=np.array([[0.0000,   0.0000,   0.0000],\
                [0.0000,   0.5000,   0.5000],\
                [0.5000,   0.0000,   0.5000],
                [0.5000,   0.5000,   0.0000]\
                ])      
    x=pos_at[:,0]
    y=pos_at[:,1]
    z=pos_at[:,2]
    return x,y,z


@mlab.show
def test_den():
    P = np.random.random((10,10))
    #P=np.load(r'c:\maxdenP.np.npy')
    fig=mlab.figure()
    
    x,y,z=gen_as()    
    pts_as=mlab.points3d(x,y,z,color=(1,0,0),colormap='gist_rainbow',figure=fig)
    x,y,z=gen_fe()  
    pts_fe=mlab.points3d(x,y,z,color=(0,1,0),colormap='gist_rainbow',figure=fig)
    x,y,z=gen_sr()  
    pts_sr=mlab.points3d(x,y,z,color=(0,0,1),colormap='gist_rainbow',figure=fig)
    im=mlab.imshow(P, colormap='gist_rainbow',figure=fig)

if __name__=='__main__':
    test_den()