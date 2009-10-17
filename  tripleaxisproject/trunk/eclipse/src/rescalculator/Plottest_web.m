% This is an example of use for ResPlot.m.
% It calculates projections of resolution ellipsoids
% in a conastant-Q scan across gapped  excitations in a 1-dimensional 
% S=1 antiferromagnet. "a" is assumed to be the chain axis. The
% dispersion is plotted using SMADemo.m.
%  ResLib v.3.3

% A. Zheludev, 1999-2006
% Oak Ridge National Laboratory

% Set up experimental conditions
setup=MakeExp; 
setup.mono.tau='PG(002)';
setup.ana.tau='PG(002)'; 
setup.mono.mosaic=30; 
setup.ana.mosaic=30; 
setup.sample.mosaic=30;
setup.sample.vmosaic=30;
setup.sample.a=6.28; 
setup.sample.b=6.28;
setup.sample.c=6.28;
setup.sample.gamma=90;
setup.hcol=[40 47 40 80];
setup.vcol=[120 120 120 120]; 
setup.orient1=[1 0 0];  
setup.orient2=[0 1 0]; 
setup.efixed=14.7
setup.infin=1;       % positive for fixed incident energy
% Set up parameters for the SMADemo.m cross section

p(1)=0;					% Gap at the AF zone-center in meV for x-axis mode
p(2)=0;					% Gap at the AF zone-center in meV for y-axis mode
p(3)=0;					% Gap at the AF zone-center in meV for z-axis mode
p(4)=30;				% Bandwidth in meV 
p(5)=0.4;				% Intrinsic Lorentzian HWHM of exccitation in meV
p(6)=1;					% Intensity prefactor
p(7)=3;                 % Flat background


% Define a constant-Q scan

%H=1.52; 
%K=0;
%L=0.0;
%W=linspace(0,20,7);

H=1;
K=0;
%K=linspace(.5,0,6);
%K=0;

%H=linspace(1,.5,6)
%K=linspace(1,1.5,6)
%H=linspace(1,.5,6)
%K=linspace(1,1.5,6)

L=0;
%H=linspace(1,.5,.1)
%K=linspace(1,1.5,.1);
%W=linspace(0,4,11);
W=[0];

%H=[H .8]
%K=[K 1.2]
%L=[L L]
%W=[W; W]

% Plot resolution parameters
setup.method=0;

EXP=setup;

[xvec,yvec,zvec,sample,rsample]=StandardSystem(EXP);
Q=modvec(H,K,L,rsample)
return
[R0C,RMC]=ResMat(Q,W,EXP)
%[R0,RMS]=ResMatS(H,K,L,W,EXP)
return


figure()
%figure(3); 
clf;
ResPlot(H,K,L,W,setup);
%ResPlotb(H,K,L,W,setup);

return

%ResPlot3D(H,K,L,W,setup,'SMADemo',p);
%ResPlot3D(H,K,L,W,setup,RANGE,'Yellow','r','r','r','SMADemo',p,SXg,SYg);

RANGE=[1.3*2*pi/setup.sample.a 1.7*2*pi/setup.sample.a 0.0*2*pi/setup.sample.a 0.7*2*pi/setup.sample.a 0 4];


% Define a mesh of (x,y) points for ploting the dispersion surface.
SX=linspace(RANGE(1),RANGE(2),100);
SY=linspace(RANGE(3),RANGE(4),40);
[SXg,SYg]=meshgrid(SX,SY);
   
% Mask one quadrant in the mesh to "cut open" the dispersion surface for a better visibility
% of resolution ellipsoids
%SXg(SXg<1.5*2*pi/3.8 & SYg<-0.15)=NaN;

ResPlot3D(H,K,L,W,setup,RANGE,'Yellow','r','r','r');

colormap cool;
