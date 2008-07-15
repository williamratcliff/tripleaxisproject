import numpy as N


def CleanArgs(**args):
    """Takes as input a set of arrays and/or scalars.  It returns a set of the arrays that are the length as the longest entry
    For shorter arrays, it pads them with the last value in the array.  It returns a dictionary based on the calling args"""
    npts=[]
    for name,val in args.iteritems():
        if type(val) in [type(13),type(13.0)]:
            args[name]=N.array([val],'float64')
        if type(val)==type([1,2,3]):
            npts.append(len(val))
        else:
            npts.append(args[name].shape[0])
    maxpts=max(npts)
    for name,val in args.iteritems():
        if type(val)==type([1,2,3]):
            for i in range(maxpts-len(val)):
                args[name].append(val[-1])
        else:
            if val.shape[0]<maxpts:
                last_val=val[-1]
                if len(val.shape)==1:
                    addendum=last_val*N.ones((maxpts-val.shape[0],),'float64')
                    args[name]=N.concatenate((val,addendum))
                elif len(val.shape)==2:
                    addendum=N.tile(last_val,(maxpts-val.shape[0],1))
                    args[name]=N.concatenate((val,addendum))
            #print name, len(val.shape),args[name]
    return args


if __name__=='__main__':
    a=N.array([1])
    b=N.array([1,2])
    c=90
    orient1=N.array([[1,2,3]])
    d=[1]
    myinput={'a':a,'b':b,'c':c}
    myoutput=CleanArgs(a=a,b=b,c=c,orient1=orient1,d=d)
    print a,b,c,orient1
    print myoutput['d']
    #print a2,b2,c2


