function y=orderparameter(p,x)
%ORDERPARAMETER  Order Parameter Function  [NEUTRON]
%
% OrderParameter
% parameters = [Amplitude
%               Tc
%               delta Tc
%               Beta
%               Background]
%
%Modified 11/26/01 by Patrick K. Mang


if (p(3)<1e-10)
    y=compute_function(p,x);
else
   % y=test_conv(p,x); for some reason using this one gives higher error (10x) in Tc
    y=cvgs(p,x);
end




%-------------------------------------
%juste compute ordinary orderparameter function
function y=compute_function(p,x)
	y=p(1)*(abs(x./p(2) - 1)).^(2*p(4));
	y=y.*(x<p(2));  %for points x>Tc set y=0
	y=y+p(5);       %add background
    
     
%------------------------------
function y=cvgs(p,x)
%convolve signal with gaussian of width p
%legacy code from Simon

xx=sort(x);   %Legacy code from Simon
%deltaxx=diff(xx);
%deltaxx=min(deltaxx(find(deltaxx~=0))); 

deltaxx=0.5;  %function requires equal step size so force step size to 0.5K

%extend limits of range to eliminate funny effects of convolution at the 
%ends of the range coming from the convolution
xx=(xx(1)-100*deltaxx):deltaxx:(xx(length(xx))+100*deltaxx);  

%this is all because of requirement that convgauss have equal step size
yy=convgauss(compute_function(p,xx),p(3)/deltaxx);
y=interp1(xx,yy,x); 

%-------------------------------------------
function z=convgauss(y,w)
%CONVGAUSS(y,w): Convolution of signal y with a gaussian of width w. It is 
%                assumed that the signal is collected at equal spacing and
%                w is in the units of the spacing.  
%Writted by Simon Larochelle

y=y(:);
lg=length(y);
nn=2^ceil(log2(lg));
y=[y ;zeros(nn-lg,1)];

p=[1 nn/2+1 w 0 0];
g=gaussian(p,1:nn)';
g=[g(nn/2+1:nn) ;g(1:nn/2)];
g=g/sum(g);

z=real(ifft(fft(g).*fft(y)));
z=z(1:lg);
    


%--------------------------------------------------------
function y_convoluted=test_conv(beta,t)
%code taken from function test_conv modified for orderparameter measurement.

dx=.5; %set grid size

hbeta=[1 100 beta(3) 0 0]; %define resolution function parameters as gaussian
x=(100-hbeta(3)*5):dx:(100+hbeta(3)*5);  %calculate h out to 5 sigma in steps of dx
if(rem(length(x),2) == 0)  %if x has an even number of elements
    x=[x x(end)+dx];        %add one more to make it an odd number
end

h=gaussian(hbeta,x);  %calculate the resolution function.

x=(t(1)-hbeta(3)*5):dx:(t(end)+hbeta(3)*5); %t defines the limit of the data we're interested in. expand the range
                           %by 5 sigma on either side to eliminate "edge effects" in the convolution.
if(rem(length(x),2) == 0)  %if x has an even number of elements
    x=[x x(end)+dx];        %add one more to make it an odd number
end

%y=compute_function(beta,x); %calculate model cross-section



z=conv(h,compute_function(beta,x));  %convolve model cross-section with resolution
n=(length(h)-1)/2; %if h is n+1 pts n:1:n then y has been extended n points on either side
z=z(n+1:length(z)-n)/sum(h); %so trim z to be the same size as y

y_convoluted=interp1(x,z,t); %interpolate to get the value of y_convoluted at actual data points t



