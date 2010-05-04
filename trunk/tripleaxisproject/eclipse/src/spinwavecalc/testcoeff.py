import sympy
import numpy as N

def coeff(expr, term):
   expr = sympy.collect(expr, term)
   symbols = list(term.atoms(sympy.Symbol))
   w = sympy.Wild("coeff", exclude=symbols)
   m = expr.match(w*term+sympy.Wild("rest"))
   if m:
       return m[w]

def exp2Trig(expr):
    

if __name__=="__main__":
    #a=sympy.Symbol('a')
    #b=sympy.Symbol('b')
    #coeff=a*b
    sympy.var("a b c d e f A B C")
    x=sympy.Symbol('x')
    y=sympy.Symbol('y')
    H = A*x**2 + b*x + c*y + d*x*y
    #expanded=expr.expand()
    #print 'expanded',expanded
    #p=sympy.Wild('p',exclude=[a,b,coeff])
    #q=sympy.Wild('q',exclude=['a','b','a*b'])
    #r=sympy.Wild('r',exclude=['a','b','a*b'])
    #pref=expanded.match(p*a**2+q*b**2+r*coeff)
    H2=a*sympy.exp(4*x)+A*sympy.exp(-4*y)
    print exp2Tri(H2)
    #print pref
    #print 'collect', sympy.collect(expanded,a)
    print coeff(H,y)
    MA = sympy.Matrix([[1,x], [y,1]])
    print MA
    print MA[0,0]
    XdX=N.zeros((2,2))
    print XdX[0,0]

    