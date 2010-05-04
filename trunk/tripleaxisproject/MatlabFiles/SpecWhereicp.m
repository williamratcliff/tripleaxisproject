function [H,K,L,E,Q,Ei,Ef]=SpecWhere(M2,S1,S2,A2,EXP)
%===================================================================================
%  function [H,K,L,E,Q,Ei,Ef]=SpecWhere(M2,S1,S2,A2,EXP)
%  ResLib v.3.3
%===================================================================================
%
%  For given values of M3,S1,S2 and A2 spectrometer motors (AKA M2,M3,M4 and M6)
%  and spectrometer and sample parameters specified in EXP calculates the wave vector
%  transfer in the sample (H, K, L), Q=|(H,K,L)|, energy tranfer E, and incident
%  and final neutron energies.
%
% A. Zheludev, 1999-2006
% Oak Ridge National Laboratory
%====================================================================================

[len,M2,S1,S2,A2,EXP]=CleanArgs(M2,S1,S2,A2,EXP);

mono=[EXP.mono];
for ind=1:len taum(ind)=GetTau(mono(ind).tau);end;
ana=[EXP.ana];
for ind=1:len taua(ind)=GetTau(ana(ind).tau); end;

ki=taum./sqrt(2-2*cos(M2));
Ei=2.072142*ki.^2;
kf=taua./sqrt(2-2*cos(A2));
Ef=2.072142*kf.^2;
E=Ei-Ef;
Q=sqrt(ki.^2+kf.^2-2*ki.*kf.*cos(S2));

Q

[orienta,orientb,z,lattice,rlattice]=StandardSystem(EXP);

%phi=-atan2(-kf.*sin(S2), ki-kf.*cos(S2)); %Angle from ki to Q
delta=abs(acos( (Q.^2+ki.^2-kf.^2)./(2*ki.*Q)));
psi=S1+delta-pi/2; %Angle from first orienting vector to to Q



qx=Q.*cos(psi);
qy=Q.*sin(psi);



H=qx*orienta(1)+qy*orientb(1);
K=qx*orienta(2)+qy*orientb(2);
L=qx*orienta(3)+qy*orientb(3);



