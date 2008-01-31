import numpy as N

def max(a,b):
    return (a<b).choose(a,b)

def min(a,b):
    return (a>b).choose(a,b)

def lookup(a,b):
    return a.searchsorted(b)

def cat(*args):
    return N.concatenate(args)

# f, f', f'', f''' = bspline3(knot, control, t, nderiv=0)
#   Evaluate the B-spline specified by the given knot sequence and
#   control values at the parametric points t.   The knot sequence
#   should be four elements longer than the control sequence.
#   Returns up to p(t), p'(t), p''(t), p'''(t) depending on nderiv.
# bspline3(knot, control, t, clamp=True)
#   Clamps the spline to the value of the final control point beyond
#   the ends of the knot sequence.   Default is 'zero' for clamping 
#   the spline to zero.
def bspline3(knot, control, t, clamp=False, nderiv=0):
  degree = len(knot) - len(control);
  if degree != 4: raise ValueError, "must have two extra knots at each end"

  if clamp:
    # Alternative approach spline is clamped to initial/final control values
    control = cat([control[0]]*(degree-1), control, [control[-1]])
  else:
    # Traditional approach: spline goes to zero at +/- infinity.
    control = cat([0]*(degree-1), control, [0])


  # Deal with values outside the range
  valid = (t > knot[0]) & (t <= knot[-1])
  tv = t[valid]
  f = N.zeros(t.shape)
  df = N.zeros(t.shape)
  d2f = N.zeros(t.shape)
  d3f = N.zeros(t.shape)
  f[t<=knot[0]] = control[0]
  f[t>=knot[-1]] = control[-1]

  # Find B-Spline parameters for the individual segments 
  end = len(knot)-1
  segment = lookup(knot,tv)-1
  tm2 = knot[max(segment-2,0)]
  tm1 = knot[max(segment-1,0)]
  tm0 = knot[max(segment-0,0)]
  tp1 = knot[min(segment+1,end)]
  tp2 = knot[min(segment+2,end)]
  tp3 = knot[min(segment+3,end)]

  P4 = control[min(segment+3,end)]
  P3 = control[min(segment+2,end)]
  P2 = control[min(segment+1,end)]
  P1 = control[min(segment+0,end)]

  # Compute second and third derivatives
  if nderiv > 1: 
    # First derivative is available almost for free
    # Second or more derivative requires extra computation
    Q4 = (P4 - P3) * 3 / (tp3-tm0)
    Q3 = (P3 - P2) * 3 / (tp2-tm1)
    Q2 = (P2 - P1) * 3 / (tp1-tm2)
    R4 = (Q4 - Q3) * 2 / (tp2-tm0)
    R3 = (Q3 - Q2) * 2 / (tp1-tm1)
    S4 = (R4 - R3) * 1 / (tp1-tm0)

    R4 = ( (tv-tm0)*R4 + (tp1-tv)*R3 ) / (tp1 - tm0)

    d2f[valid] = R4
    d3f[valid] = S4

  # Compute function value and first derivative
  P4 = ( (tv-tm0)*P4 + (tp3-tv)*P3 ) / (tp3 - tm0)
  P3 = ( (tv-tm1)*P3 + (tp2-tv)*P2 ) / (tp2 - tm1)
  P2 = ( (tv-tm2)*P2 + (tp1-tv)*P1 ) / (tp1 - tm2)
  P4 = ( (tv-tm0)*P4 + (tp2-tv)*P3 ) / (tp2 - tm0)
  P3 = ( (tv-tm1)*P3 + (tp1-tv)*P2 ) / (tp1 - tm1)
  fastdf = (P4-P3) * 3 / (tp1-tm0)
  P4 = ( (tv-tm0)*P4 + (tp1-tv)*P3 ) / (tp1 - tm0)

  # Check that fast df calculation matches the direct Q4 calculation.
  # if nderiv > 1: print "|fast df - df| = ",norm(df-Q4)

  df[valid] = fastdf
  f[valid] = P4

  if nderiv == 0: return f
  elif nderiv == 1: return f,df
  elif nderiv == 2: return f,df,d2f
  else: return f,df,d2f,d3f

# Assertions left over from original octave code --- I'm not ready
# to write a generic assert yet in Python
#!assert(bspline3([0 0 0 1 1 3 4 6 6 6],[0 0 0 0 0 0],2.2),0,10*eps);
#!assert(bspline3([0 0 0 1 1 3 4 6 6 6],[1 1 1 1 1 1],2.2),1,10*eps);
#!assert(bspline3([0 0 0 0 1 4 5 5 5 5],[1:6],2),761/240,10*eps);
#!assert(bspline3([0 0 0 0 1 4 5 5 5 5],[1:6],[2,2]),[761/240,761/240],10*eps);
#!assert(bspline3([0 0 0 1 1 3 4 6 6 6],[1:6],3.2),4.2976,10*eps);

import numpy as nx

class BSpline3:
    """Manage control points for parametric B-spline."""
    # TODO: this class doesn't give much control over knots.
    def __init__(self, x, y, clamp=True):
        n = len(x)
        self.knot = nx.concatenate([[0.]*2, range(n), [n-1]*2])
        self.x = x
        self.y = y
        self.clamp = clamp
        
    def __len__(self): 
        """Count the knots"""
        return len(self.x)

    def __getitem__(self, i):
        """Set control point for a knot"""
        return self.x[i], self.y[i]
    
    def __setitem__(self, i, pair):
        """Get control point for a knot"""
        self.x[i],self.y[i] = pair
    
    def __delitem__(self, i):
        """Delete a knot"""
        if i < 0 or i >= len(self.x): raise IndexError
        self.x = nx.delete(self.x,i)
        self.y = nx.delete(self.y,i)
        self.knot = nx.delete(self.knot,i+2)
        if i == 0:
            self.knot[0:2] = self.knot[2]
        elif i == len(self.x)-2:
            self.knot[-2:-1] = self.knot[-3]
            
    def __call__(self, t):
        """Evalaute a B-spline at points t"""
        fx = bspline3(self.knot,self.x,t,clamp=self.clamp)
        fy = bspline3(self.knot,self.y,t,clamp=self.clamp)
        return fx,fy

    def append(self,x,y):
        """Add a knot to the end"""
        self.x = nx.concatenate([self.x,[x]])
        self.y = nx.concatenate([self.y,[y]])
        k = self.knot[-1]+1
        self.knot = nx.concatenate([self.knot,[k]])
        self.knot[-3:-1] = k
        
    def sample(self,n=400):
        """Sample the B-spline at n equidistance points in t"""
        t = nx.linspace(self.knot[2],self.knot[-3],n)
        return self.__call__(t)
        

def demo():
  import pylab
  t = N.linspace(-1,7,40 );
  knot = N.array([0, 1, 1, 3, 4, 6],'f')
  #knot = N.array([0, 0, 1, 4, 5, 5],'f')
  control = N.array([1, 2, 3, 2, 1, 2],'f') 
  knotseq = cat([knot[0]-1,knot[0]], knot, [knot[-1],knot[-1]+1])
  f = bspline3(knotseq,control,t,clamp=True);
  #print zip(t,f)
  pylab.plot(t,f,'-',knot,control,'x');
  pylab.show()

if __name__ == "__main__": demo()
