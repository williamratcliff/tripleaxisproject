import sympy
I=sympy.I

if __name__=="__main__":
    x=sympy.Symbol('x',real=True)
    expr=sympy.exp(I*x)
    print expr.expand(complex=True)