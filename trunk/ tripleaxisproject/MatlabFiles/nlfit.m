function [beta,chisquare,errors,fitresult] = nlfit(X,y,model,beta0,varargin)

%NLFIT Nonlinear least-squares data fitting by the Levenberg-Marquardt method.
%   NLFIT(X,Y,'MODEL',BETA0,varargin) finds the coefficients of the nonlinear 
%   function described in MODEL. MODEL is a user supplied function having 
%   the form y = f(beta,x). That is MODEL returns the predicted values of y
%   given initial parameter estimates, beta, and the independent variable, X.   
%
%   [BETA,CHISQUARE,ERRORS,FITRESULT] = NLINFIT(X,Y,'MODEL',BETA0,varargin) 
%   returns the fitted coefficients BETA, the goodness-of-fit measure chisquare,
%   the errors on each fitted parameter in BETA, and a structure FITRESULT, 
%   containing more detailed information on the fit.
%
%   The arguments varargin are optional. The first possible arguement is a vector 
%   specifing which parameters are to be fit via ones and zeros. If an empty matrix 
%   is passed or no arguement is specified all parameters will be fit.
%
%   The second optional arguement sigma controls the weighting of the fit. It should 
%   consist of a vector the size of Y containing the errorbar on each Y point. If 
%   this arguement is not specified by default sigma=sqrt(Y), as is appropriate for 
%   fitting scattering datasets.
%  
%   The third optional argument may be a single item or structure that is passed unchanged
%   to the model function.  This is useful if more complex additional information is necessary
%   in the function.
%
%   B.A. Jones 12-06-94.
%   Copyright (c) 1993-97 by The MathWorks, Inc.
%   $Revision: 2.9 $  $Date: 1997/04/08 15:09:39 $
%
%
%Modified 3/04/99 by Simon Larochelle
%Modified 2/01/01 by Patrick K. Mang
%Modified 11/19/01 by Simon Larochelle


IterationPrint=1;

%initialize optional input variables. a is defined so as to fit all parameters.
%sigma, the weighting in the fit is defined for statistical weighting.
a=ones(size(beta0));
sig=sqrt(y(:));
%sig=ones(size(y));

%now if any optional parameters are specified, set a and sigma to the 
%appropriate values
switch nargin  
case 5
   a=deal(varargin{1});
case 6
   [a,sig]=deal(varargin{:});
   if isempty(a)
      a=ones(size(beta0)); %default to fitting all parameters
   end
case 7
    [a,sig,fun_struct]=deal(varargin{:}); %any additional parameters are passed on to model function 
    if isempty(a)                   %as a structure
        a=ones(size(beta0));
        sig=sqrt(y(:));
    end
end

 
n = length(y);  %define number of pts
if min(size(y)) ~= 1
   error('Requires a vector second input argument.');
end
X = X(:);  %columnize vectors
y = y(:);
beta0 = beta0(:);  
sig=sig(:);
sig(find(sig==0))=inf;  %guard against division by zero

A=find(a(:));  %A defines which parameters are being varied. eg A=[1 2 4];
p=length(A);   %define number of parameters being fitted.

J = zeros(n,p);    %initialize Jacobian
beta = beta0;      %initialize parameters
betanew=beta0;     
betaLM=beta0;
betanew(A) = beta(A)*1.01; %slightly vary the parameters.  A serves to ...
                          %pick out only those parameters being varied.

                          
if strcmp(model,'twoaxis')
    maxiter=10;
else
    maxiter = 100;
end

iter = 0;
betatol = 1e-4; %fractional tolerance of variations in parameters
rtol = 1e-4; %fractional tolerance on chi-square
sse = 1;     %sse is Chi-square
sseold = 1;
lambda=0.01; 


while (~isempty(find((abs((betanew-beta)./(beta+sqrt(eps))) > betatol))) | ...
      ((sseold-sse)/(sse+sqrt(eps)) > rtol)) & ...
      (iter < maxiter)
   
   %while either the fractional change in any parameter beta or chi-square is 
   %greater than specified tolerances AND the number of iterations 
   %is less than max
    
    if iter > 0,
       beta = betanew;
    end
    
iter = iter + 1;

if IterationPrint
    fprintf('Iteration #%d \t Chi-Square %f\n',iter,sse/(n-p));
end
   
if nargin > 6
   yfit = feval(model,beta,X,fun_struct);   %evaluate the model
else
    yfit=feval(model,beta,X);
end

   r = (y - yfit)./sig;          %calculate residual
   sseold = r'*r;                %calculate chi-square

   for k = 1:p
      %calculate J=dy/da
      if nargin > 6
          J(:,k)=nlfit_deriv(model,beta,X,yfit,A(k),nargin,fun_struct)./sig;
      else
          J(:,k)=nlfit_deriv(model,beta,X,yfit,A(k),nargin)./sig;
      end
      %A(k) cleverly passes wrt which parameter in beta should the derivative
      %be calculated as k cycles through the total number of varied parameters. 
      %This step calculates the derivative of the model, yfit, wrt to small changes
      %in the parameters, beta. It calculates it at each point, returning #pts 
      %derivatives.
   end


   % Levenberg-Marquardt type adjustment 
   % Gauss-Newton step -> J\r
   % LM step -> inv(J'*J+constant*eye(p))*J'*r

   JJ=J'*J;  %calculate J squared which is equal to alpha_ll
   Jr=J'*r;  %calculate dChi-square/dbeta(n) = beta						
   stepLM = (JJ-diag(diag(JJ))+JJ.*(eye(p)*(1+lambda)))\Jr;
   %first step replaces the diagonal elements of JJ, which are equal to alpha_ll
   %with (1+lambda)*alpha_ll
   %Dividing by Jr, which is equal to beta, then computes the next step of the parameters
   %according to LM algorithim. c.f. 15.5.14 in Numerical Recipes in C.
   
   betaLM(A) = beta(A) + stepLM;
   if nargin>6
       yfitnew = feval(model,betaLM,X,fun_struct);
   else
       yfitnew = feval(model,betaLM,X);
   end
   rnew = (y - yfitnew)./sig;
   sseLM = rnew'*rnew;    %calculate new chi-square after moving by LM step
   iter1 = 0;
   while sseLM > sseold & iter1 < 12   %keep stepping until we reduce chi-square
      stepLM = stepLM/sqrt(10);
      betaLM(A) = beta(A) + stepLM;
      if nargin > 6
          yfitnew = feval(model,betaLM,X,fun_struct);
      else
          yfitnew = feval(model,betaLM,X);
      end
      rnew = (y - yfitnew)./sig;
      sseLM = rnew'*rnew;
      iter1 = iter1 + 1;
   end

    if iter1 < 12
      lambda=lambda/2;
      betanew=betaLM;
      sse=sseLM;
	else 
      lambda=lambda*10.;
   end
 end   


   


if iter == maxiter
   disp('NLINFIT did NOT converge. Returning results from last iteration.');
end

%main fitting routine is over. now prepare output values
chisquare=sse/(n-p);   %normalize chi-square by degrees of freedom
fitresult.iterations = iter;
fitresult.phi = sse;
fitresult.lambda = lambda;  %prepare an optional structure fitresult
fitresult.NumPts = n;
fitresult.NumParam = p;
fitresult.yfit = yfit;
fitresult.residual = r;
fitresult.sigma = sig;
fitresult.J=J;
errors=nlparci(beta(A),r,J,A); %calculate errors in fit only for parameters being fit







function y=nlfit_deriv(model,beta,X,Y,n,nargin_nlfit,fun_struct)
%nlfit_deriv(beta,x,n): Take the first derivative of y with respect to beta(n) at
%			points X. This is a subroutine of nlfit.
%        
%        nlfit_deriv only calculates dy/dbeta at each point. Later in the main 
%        routine they will all have to be summed and weighted together.

delta = zeros(size(beta));  %initialize delta to zero
delta(n) = sqrt(eps)*beta(n); %make delta for the parameter we are calculating very
if delta(n)==0                %small. All others are still zero.
        delta(n)=eps;
end

if nargin_nlfit > 6
    y1 = feval(model,beta+delta,X,fun_struct);
else
    y1 = feval(model,beta+delta,X);
end

% y2 = feval(model,beta,X);
y = (y1 - Y)/delta(n);  %calculate dy/dbeta
%y = (y1 - Y)/(sqrt(eps)*beta(n));  

while(sum(y)==0 & delta(n)<0.01*beta)
        delta(n)=delta(n)*10.;
        if nargin_nlfit > 6
            y1=feval(model,beta+delta,X,fun_struct);
        else
            y1=feval(model,beta+delta,X);
        end
        y=(y1-Y)/(sqrt(eps)*beta(n));
end

%end of nlfit_deriv




function delta = nlparci(x,f,J,A)
%code with slight modifications from function nlparci.m
%returns one gaussian standard deviation (confidence interval
%equals 68%) estimate of errors on fitted parameters.

%NLPARCI Confidence intervals on parameters of nonlinear models.
%   CI = NLPARCI(X,F,J) returns the 95% confidence interval CI
%   on the nonlinear least squares parameter estimate X, given the 
%   residual sum of squares, F, and the Jacobian matrix ,J, at the solution.
%
%   The confidence interval calculation is valid for systems where 
%   the number of rows of J exceeds the length of X. 
%
%   NLPARCI uses the outputs of NLINFIT for its inputs.
%   Example:
%      [x,f,J]=nlinfit(input,output,model,xinit);
%      ci = nlparci(x,f,J);
%
%   See also NLINFIT.
%

%   Bradley Jones 1-28-94
%   Copyright (c) 1993-98 by The MathWorks, Inc.
%   $Revision: 2.7 $  $Date: 1998/07/10 14:45:54 $

%initialization
if nargin < 3
   error('Requires three inputs.');
end;

f = f(:);
[m,n] = size(J);
if m <= n
   error('The number of observations must exceed the number of parameters.');
end;

if length(x) ~= n
   error('The length of x must equal the number of columns in J.')
end

% approximation when a column is zero vector
temp = find(max(abs(J)) == 0);
if ~isempty(temp)
   J(temp,:) = J(temp,:) + sqrt(eps);
end;

%calculate covariance
[Q R] = qr(J,0);
Rinv = R\eye(size(R));
diag_info = sum((Rinv.*Rinv)')';

v = m-n;
rmse = sqrt(sum(f.*f)/v);

% calculate one standard deviation
%note that the number .8413 is .5(1+.68), the proper
%formula for calculating a 68% confidence integral
%using the t-distribution.

delta(A) = sqrt(diag_info) .* rmse*tinv(0.8413,v); 
%ci = [(x(:) - delta) (x(:) + delta)];

%--end of nlparci.m---





