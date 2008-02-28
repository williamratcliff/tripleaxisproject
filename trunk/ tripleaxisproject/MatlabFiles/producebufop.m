function producebufop()
Ei=35
% scattering plane has to be a-b plane (redefine lattice directions if necessary)
% (h,2h,l) scattering zone

lowerA3=-80;upperA3=100;
Qord1=[0.5;+1.0;0]; % first entry is (m,0,0), second is (0,0,1)
% Qord1=[0.5;0;0];
u=1;
ki=sqrt(Ei/2.072);kf=sqrt(Ei/2.072);
% deff=sqrt(1/((1/7.294)^2+(2/8.551)^2));
astar=2*pi/5.6719;bstar=2*pi/5.6818;cstar=1;
for p2=1:2:1
    for m=-11:1:11
        for n=-11:1:11
            Q=[m;n;0]+p2*Qord1;  % m specified (m,2*m,0); n specifies (0,0,n)
             % [m,n,0] == (m,2*m,n)
            if norm(Q)~=0
                [psi theta beta alpha]=hkltoa3a4(Q,ki,kf,astar,bstar,cstar);
                A3=90-psi/pi*180;A4=theta/pi*180;
                if A3>lowerA3 & A3<upperA3 & (A4>5) & (A4<110)
                    A3m(u)=A3;A4m(u)=A4;ss(u)=round(A4/30);tt(u)=round(A3/40);
                    % x(u)=m+p2*0.5;y(u)=2*(m+p2*0.5);z(u)=n+p2*0.25;
                    x(u)=m+p2*Qord1(1);y(u)=0;z(u)=n+p2*Qord1(2);
%                     if mod(n+m,2)==0, mon(u)=20000;,elseif mod(n+m+1,2)==0,mon(u)=40000;,end
                    mon(u)=30000;
                    u=u+1;
                end
            end
        end
    end
end
format short;
Res=[A3m' A4m' ss' tt' x' y' z' mon'];Ressort=Res;
[Ressort,I]=sortrows(Res,[4 3]);
% [Ressort,I]=sortrows(Ressort,[5 6 7]);
for v=1:1:u-1
    fprintf(1,'BUFOP b3 hc=%3.3f, kc=%2.3f, lc=%2.3f \nrb3 \n',Ressort(v,5),Ressort(v,6),Ressort(v,7));%, MPF=30
    fprintf(1,'BUFOP I1 A3=%3.3f, A4=%2.3f \n\n',Ressort(v,1),Ressort(v,2));
%    fprintf(1,'BUFOP I1 A3=%3.3f, A4=%2.3f, I3=0.07, pts=51\nri1 \n',Ressort(v,1),Ressort(v,2));
end
disp('  min(A3m)  max(A3m)  min(A4m)  max(A4m)');
disp([min(A3m) max(A3m) min(A4m) max(A4m)]);
[numRef, dumb]=size(Res);fprintf(1,'  Number of accessible reflections: %3.0f \n',round(numRef));