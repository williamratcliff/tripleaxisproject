function [Ifx,errbnd] = quadva(f,interval,reltol,abstol)
% When INTERVAL = [A,B], computes the integral of a continuous function 
% f(x) from A to B.  A can be -Inf and/or B can be Inf. F accepts a row
% vector X and returns a vector Y with Y(m) = f(X(m)) for m = 1,... .
%
% QUADVA returns an approximation Ifx to the integral and optionally an 
% approximate bound, ERRBND, on the error |integral - Ifx|. It attempts 
% to compute Ifx such that |Ifx - integral| <= max(ABSTOL,RELTOL*|Ifx|).
% If QUADVA is unsuccessful, Ifx and ERRBND are still meaningful, so a 
% warning is issued that includes ERRBND.   
%
% If the interval is infinite, say [a,Inf), then for the integral to 
% exist, f(x) must decay as x -> Inf and QUADVA requires it to decay 
% rapidly. Special methods should be used for oscillatory functions on 
% infinite intervals, but QUADVA can be used if f(x) decays fast enough.
%
% QUADVA will integrate f(x) that are singular at finite end points if
% the singularities are not too strong. For example, it will integrate
% f(x) that behave like log|x-c| or |x-c|^p for p >= -1/2 with c = a 
% and/or c = b. If f(x) is singular at points inside (A,B), write the 
% integral as a sum of integrals over subintervals with the singular
% points as end points, compute them with QUADVA, and add the results.
%
% QUADVA starts with samples of f(x) at 150 points in (A,B). It must be
% able to recognize the behavior of f(x) from these samples, so if f(x) 
% oscillates very rapidly or has sharp peaks, it may be necessary to 
% subdivide the interval. To do this, make INTERVAL an array with entries
% (breakpoints) that increase from A = INTERVAL(1) to B = INTERVAL(end). 
% If f(x) is only piecewise smooth, the points of discontinuity should
% be breakpoints. If f(x) has very sharp peaks, they should occur at
% midpoints of subintervals. 
%
% Written by L.F. Shampine, Mathematics Department, Southern Methodist
% University, lshampin@smu.edu.  All rights reserved to the author.

interval = interval(:)';
if any(diff(interval) <= 0)
    error('Entries of INTERVAL must increase.')
end
a = interval(1); b = interval(end);

% Generally the error test is a mixed one, but pure absolute error
% and pure relative error are allowed.  If a pure relative error
% test is specified, the tolerance must be at least 100*EPS. Defaults
% are 1e-5 for relative error, 1e-10 for absolute error.
if nargin < 3 || isempty(reltol)
    rtol = 1e-5;
elseif reltol <= 0
    rtol = 0;
else
    rtol = max(reltol,100*eps);
end
if nargin < 4 || isempty(abstol)
    atol = 1e-10;
else
    atol = max(abstol,0);
end
if atol+rtol == 0
    rtol = 1e-5;
    atol = 1e-10;
end

% Gauss-Kronrod (7,15) pair. Use symmetry in defining nodes and weights.
pnodes = [ 0.2077849550078985; 0.4058451513773972; 0.5860872354676911; ...
           0.7415311855993944; 0.8648644233597691; 0.9491079123427585; ...
           0.9914553711208126 ];
nodes = [ -pnodes(end:-1:1); 0; pnodes ];
samples = 15;
pwt = [ 0.2044329400752989; 0.1903505780647854; 0.1690047266392679; ...
        0.1406532597155259; 0.1047900103222502; 0.06309209262997855; ...
        0.02293532201052922 ];
wt = diag([ pwt(end:-1:1); 0.2094821410847278; pwt]);
pwt7 = [ 0; 0.3818300505051189; 0; 0.2797053914892767; 0; ...
         0.1294849661688697; 0 ];
ewt = wt - diag([ pwt7(end:-1:1); 0.4179591836734694; pwt7 ]);

warns = warning('query','all'); 
% Identify the task. If breakpoints are specified, work out 
% how they map into the standard interval.
if isfinite(a) && isfinite(b)
    task = 1;
    if length(interval) > 2
        tinterval = invert1(interval);
    else
        tinterval = linspace(-1,1,11);
    end
    [Ifx,errbnd,OK] = Vadapt(@f1,tinterval,rtol,atol);
elseif isfinite(a) && isinf(b)
    task = 2;
    if length(interval) > 2
        alpha = sqrt(interval(2:end-1) - a);
        tinterval = [0, alpha./(1 + alpha), 1];
    else
        tinterval = linspace(0,1,11);
    end
    [Ifx,errbnd,OK] = Vadapt(@f2,tinterval,rtol,atol);
elseif isinf(a) && isfinite(b)
    task = 3;
    if length(interval) > 2
        alpha = sqrt(b - interval(2:end-1));
        tinterval = [ -1, -alpha./(1 + alpha), 0];
    else
        tinterval = linspace(-1,0,11);
    end    
    [Ifx,errbnd,OK] = Vadapt(@f3,tinterval,rtol,atol);
elseif isinf(a) && isinf(b)
    task = 4;
    if length(interval) > 2
        tinterval = invert4(interval);
    else
        tinterval = linspace(-1,1,11);
    end    
    [Ifx,errbnd,OK] = Vadapt(@f4,tinterval,rtol,atol);
end
warning(warns);
if ~OK
   str = sprintf(['\n\n***Ifx does not satisfy error test. ',...
                 'Approximate bound on error is %9.1e.\n\n'],errbnd);
   fprintf(str)
end 

%==Nested functions========================================================
function [Ifx,errbnd,OK] = Vadapt(f,tinterval,rtol,atol)

% Initialize array of subintervals of [a,b].
nsubs = length(tinterval) - 1;
subs = [tinterval(1:end-1);tinterval(2:end)];

% Initialize partial sums.
IfxOK = 0;
errOK = 0;

% Initialize main loop
OK = true;
first = true;
Ifx = []; errbnd = [];
while true
    % SUBS contains subintervals of [a,b] where the integral is not
    % sufficiently accurate.  The first row of SUBS holds the left end
    % points and the second row, corresponding right end points.
    midpt = sum(subs)/2;   % midpoints of the subintervals
    halfh = diff(subs)/2;  % half the lengths of the subintervals
    x = nodes*halfh + repmat(midpt,samples,1);  
    x = reshape(x,1,[]);   % function f expects a row vector
    if ~first, warning('off','all'); end
    [fx,too_close] = feval(f,x); 
    % Quit if mesh points are too close or too close to a
    % singular point or got into trouble on first evaluation.
    if too_close || any(~isfinite(fx)), break; end
    fx = reshape(fx,samples,[]);

    % Quantities for subintervals.
    Ifxsubs = sum(wt*fx).*halfh;  
    errsubs = sum(ewt*fx).*halfh;
    
    % Quantities for all of [a,b].
    Ifx = sum(Ifxsubs) + IfxOK;      
    errbnd = abs( sum(errsubs) + errOK );

    % Test for convergence:
    if ( errbnd <= max(atol,rtol*abs(Ifx)) ), return; end
    
    % Locate subintervals where the integrals are sufficiently
    % accurate.  Update partial sums and remove these subintervals 
    % from SUBS.  Split the remaining subintervals in half.
    ndx = find(abs(errsubs) <= (2/(b-a))*halfh ...
                               .*max(atol,rtol*abs(Ifxsubs)));
    errOK = errOK + sum(errsubs(ndx));
    IfxOK = IfxOK + sum(Ifxsubs(ndx));
    subs(:,ndx) = [];  
    % Quit if splitting results in too many subintervals.
    if (2*size(subs,2) > 650), break; end
    midpt = sum(subs)/2;
    subs = reshape([subs(1,:); midpt; midpt; subs(2,:)],2,[]);
    first = false;
end
if first
    warning(warns);
    error('Difficulty evaluating integrand.');
end
OK = false;

end %Vadapt

function [y,too_close] = f1(t)
% Transform to weaken singularities at both ends: [a,b] -> [-1,1]
    Tt = 0.25*(b-a)*t.*(3 - t.^2) + 0.5*(b+a);
    too_close = check_spacing(Tt);
    if too_close
        y = [];
    else
        y = feval(f,Tt);
        y = 0.75*(b-a)*y.*(1 - t.^2);
    end
end % f1

function [y,too_close] = f2(t)
% Transform to weaken singularity at left end: [a,Inf) -> [0,Inf). 
% Then transform to finite interval: [0,Inf) -> [0,1].
    Tt = t ./ (1 - t);
    T2t = a + Tt.^2;
    too_close = check_spacing(T2t);
    if too_close
        y = [];
    else
        y = feval(f,T2t);
        y =  2*Tt .* y ./ (1 - t).^2;  
    end
end % f2

function [y,too_close] = f3(t)
% Transform to weaken singularity at right end: (-Inf,b] -> (-Inf,b]. 
% Then transform to finite interval: (-Inf,b] -> (-1,0].
    Tt = t ./ (1 + t);
    T2t = b - Tt.^2;
    too_close = check_spacing(T2t);
    if too_close
        y = [];
    else
        y = feval(f,T2t);
        y = -2*Tt .* y ./ (1 + t).^2;
    end
end % f3 

function [y,too_close] = f4(t)
% Transform to finite interval: (-Inf,Inf) -> (-1,1).        
    Tt = t ./ (1 - t.^2);
    too_close = check_spacing(Tt);
    if too_close
        y = [];
    else
        y = feval(f,Tt);
        y = y .* (1 + t.^2) ./ (1 - t.^2).^2;
    end
end % f4

function too_close = check_spacing(x)
    ax = abs(x);
    too_close = any(diff(x) <= 100*eps*max(ax(1:end-1),ax(2:end)));
end % too_close

function tint = invert1(xint)
% Invert transformation from [a,b] to [-1,+1].
% Assumes xint is a row vector and returns a row vector.
    alpha = (xint(2:end-1) - 0.5*(b+a) )/ (0.25*(b-a));
    % Solve P(t) = t^3 -3t + alpha = 0.  Note -2 < alpha < +2.
    % Converges monotonely and quadratically from guess of zero.
    % First iterate is immediately available:
    iterate = alpha/3; 
    % Newton's method: del = - P(iterate)/P'(iterate)
    % P'(t) > 0, so method is quadratically convergent.
    while true
        del = - (iterate.^3 - 3*iterate + alpha)./(3*iterate.^2 - 3);
        iterate = iterate + del;
        if all(abs(del) <= 10*eps), break; end;
    end
    tint = [-1, iterate, +1];
end % invert1

function tint = invert4(xint)
% Invert transformation from (-Inf,+Inf) to (-1,+1).
% Assumes xint is a row vector and returns a row vector.
    ndx = find(xint(2:end-1) < 0) + 1;
    if isempty(ndx)
        alphan = [];
    else
        alphan = (-1 + sqrt(1 + 4*xint(ndx).^2)) ./ ...
                 (2*xint(ndx));
    end        
    if any(xint(2:end-1) == 0)
        alpha0 = 0;
    else
        alpha0 = [];
    end
    ndx = find(xint(2:end-1) > 0) + 1;
    if isempty(ndx)
        alphap = [];
    else
        alphap = (-1 + sqrt(1 + 4*xint(ndx).^2)) ./ ...
                 (2*xint(ndx));
    end
    tint = [ -1, alphan, alpha0, alphap, +1];
end % invert4


%==========================================================================

end % quadva