import numpy as N
import lattice_calculator
import rescalc
import prefdemo
#import sqwdemo
import sqwdemo_bifeo3
import smademo
import convres
import convres_sma

if __name__=="__main__":
    if 1:
        a=N.array([6],'d')
        b=N.array([7],'d')
        c=N.array([8],'d')
        alpha=N.array([90],'d')
        beta=N.array([90],'d')
        gamma=N.array([90],'d')
 #       orient1=N.array([[0,1,1]],'d')
        orient1=N.array([[1,0,0]],'d')
        orient2=N.array([[0,1,0]],'d')
        mylattice=lattice_calculator.lattice(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,\
                               orient1=orient1,orient2=orient2)
        H=N.array([1.5,1.5],'d');K=N.array([0,0],'d');L=N.array([.35,0.35],'d');W=N.array([20,10],'d')
        EXP={}
        EXP['ana']={}
        EXP['ana']['tau']='pg(002)'
        EXP['mono']={}
        EXP['mono']['tau']='pg(002)';
        EXP['ana']['mosaic']=25
        EXP['mono']['mosaic']=25
        EXP['sample']={}
        EXP['sample']['mosaic']=25
        EXP['sample']['vmosaic']=25
        EXP['hcol']=N.array([80, 40, 40, 80],'d')
        EXP['vcol']=N.array([120, 120, 120, 120],'d')
        EXP['infix']=-1 #positive for fixed incident energy
        EXP['efixed']=14.7
        EXP['method']=0
        setup=[EXP,EXP]
        newinput=lattice_calculator.CleanArgs(a=a,b=b,c=c,alpha=alpha,beta=beta,gamma=gamma,orient1=orient1,orient2=orient2,\
                            H=H,K=K,L=L,W=W,setup=setup)
        mylattice=lattice_calculator.lattice(a=newinput['a'],b=newinput['b'],c=newinput['c'],alpha=newinput['alpha'],\
                        beta=newinput['beta'],gamma=newinput['gamma'],orient1=newinput['orient1'],\
                        orient2=newinput['orient2'])
        #myrescal.__init__(mylattice)

        #Parameter values for the cross section
        p=N.zeros((7,),'float64')
        p[0]=3#					% Gap at the AF zone-center in meV for x-axis mode
        p[1]=3#				% Gap at the AF zone-center in meV for y-axis mode
        p[2]=3#					% Gap at the AF zone-center in meV for z-axis mode
        p[3]=30#				% Bandwidth in meV
        p[4]=0.4#				% Intrinsic Lorentzian HWHM of exccitation in meV
        p[5]=6e4#					% Intensity prefactor
        p[6]=40#                 % Flat background
        #print 'p ',p
        myrescal=rescalc.rescalculator(mylattice)
        (prefactor,background)=prefdemo.PrefDemo(H,K,L,W,myrescal,p)
        #print 'prefactor ',prefactor
        #print 'background ',background
        sqw=sqwdemo_bifeo3.SqwDemo(H,K,L,W,p)
        disp,myint,WL=smademo.SMADemo(H,K,L,p)
        #print disp
        #print myint
        #print WL
        #print 'sqw ', sqw
        ac=[5,0]
        conv=convres.ConvRes(sqwdemo_bifeo3.SqwDemo,prefdemo.PrefDemo,H,K,L,W,myrescal,setup,p,METHOD='fixed',ACCURACY=ac)
        print conv
        exit()
        conv_sma=convres_sma.ConvResSMA(smademo.SMADemo,prefdemo.PrefDemo,H,K,L,W,myrescal,setup,p,METHOD='fixed',ACCURACY=ac)
        print conv_sma
        #R0,RMS=myrescal.ResMatS(H,K,L,W,setup)
        #myrescal.ResPlot(H, K, L, W, setup)
        #print 'RMS'
        #print RMS.transpose()[0]
        #print myrescal.calc_correction(H,K,L,W,setup,qscan=[[1,1,0],[1,1,0]])
        #print myrescal.CalcWidths(H,K,L,W,setup)