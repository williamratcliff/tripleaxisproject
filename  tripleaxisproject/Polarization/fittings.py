from scipy import optimize
#from numpy import *



class Parameter:
    """
        Class to handle model parameters
    """
    def __init__(self, model, name, value=None):
            self.model = model
            self.name = name
            if not value==None:
                self.model.setParam(self.name, value)
           
    def set(self, value):
        """
            Set the value of the parameter
        """
        self.model.setParam(self.name, value)

    def __call__(self):
        """ 
            Return the current value of the parameter
        """
        return self.model.getParam(self.name)
    
def sansfit(model, pars, x, y, err_y ,qmin=None, qmax=None):
    """
        Fit function
        @param model: sans model object
        @param pars: list of parameters
        @param x: vector of x data
        @param y: vector of y data
        @param err_y: vector of y errors 
    """
    def f(params):
        """
            Calculates the vector of residuals for each point 
            in y for a given set of input parameters.
            @param params: list of parameter values
            @return: vector of residuals
        """
        i = 0
        for p in pars:
            p.set(params[i])
            i += 1
        
        residuals = []
        for j in range(len(x)):
            if x[j]>qmin and x[j]<qmax:
                residuals.append( ( y[j] - model.runXY(x[j]) ) / err_y[j] )
       
        return residuals
        
    def chi2(params):
        """
            Calculates chi^2
            @param params: list of parameter values
            @return: chi^2
        """
        sum = 0
        res = f(params)
        for item in res:
            sum += item*item
        return sum
        
    p = [param() for param in pars]
    out, cov_x, info, mesg, success = optimize.leastsq(f, p, full_output=1, warning=True)
    print info, mesg, success
    # Calculate chi squared
    if len(pars)>1:
        chisqr = chi2(out)
    elif len(pars)==1:
        chisqr = chi2([out])
        
    return chisqr, out, cov_x    


def calcCommandline(self,event):
    """
        Testing implementation
    """
  
    # Fit a Line model
    from LineModel import Line
    line    = Line()
    cstA = Parameter(line, 'A', event.cstA)
    cstB  = Parameter(line, 'B', event.cstB)        
    y = line.run()
    chisqr, out, cov = sansfit(line, [cstA, cstB],  event.x, y, 0) 
    # print "Output parameters:", out
    print "The right answer is [70.0, 1.0]"
    print chisqr, out, cov

if __name__ == "__main__": 
    
    calcCommandline()
