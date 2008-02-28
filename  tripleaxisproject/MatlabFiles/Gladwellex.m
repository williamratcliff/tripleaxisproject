function Gladwellex
% Examples from I. Gladwell's MA3315, Fall 2005 class.

reltol = 0;
abstol = 1e-12;

disp('============================================================')
disp(' ')
disp('Examples from I. Gladwell''s MA3315, Fall 2005 class.')
fprintf('Absolute error tolerance of %g.\n',abstol);
disp(' ')
disp('============================================================')
disp(' ');
disp('Example 1')
Itrue = -1/3*exp(-24)+1/3; % Using symbolic capabilities
disp(' ')
disp('Using QUADVA')
nvfe = 0;
If = quadva(@f1,[0,8],reltol,abstol);
fprintf('Approximate value of the integral = %f\n',If);
fprintf('True error        %e\n',abs(Itrue - If));
fprintf('%i (vector) evaluations of integrand were used.\n\n',nvfe);
disp('Using QUADL')
nvfe = 0;
If = quadl(@f1,0,8,abstol);
fprintf('Approximate value of the integral = %f\n',If);
fprintf('True error is      %e\n',abs(Itrue - If));
fprintf('%i (vector) evaluations of integrand were used.\n\n',nvfe);

disp('============================================================')
disp(' ');
disp('Example 2')
Itrue = 35/6+1/2*2^(1/2)-1/3*3^(1/2); % Using symbolic capabilities
disp(' ')
disp('Using QUADVA without breakpoints')
nvfe = 0;
If = quadva(@f2,[-1,2],reltol,abstol);
fprintf('Approximate value of the integral = %f\n',If);
fprintf('True error        %e\n',abs(Itrue - If));
fprintf('%i (vector) evaluations of integrand were used.\n\n',nvfe);
disp('Using QUADL')
nvfe = 0;
If = quadl(@f2,-1,2,abstol);
fprintf('Approximate value of the integral = %f\n',If);
fprintf('True error is      %e\n',abs(Itrue - If));
fprintf('%i (vector) evaluations of integrand were used.\n\n',nvfe);

disp('Using QUADVA with breakpoints')
nvfe = 0;
If = quadva(@f2,[-1,-1/sqrt(2),1/sqrt(3),2],reltol,abstol);
fprintf('Approximate value of the integral = %f\n',If);
fprintf('True error        %e\n',abs(Itrue - If));
fprintf('%i (vector) evaluations of integrand were used.\n\n',nvfe);

disp('============================================================')
disp(' ');
disp('Example 3')
Itrue = 3;
disp(' ')
disp('Using QUADVA')
nvfe = 0;
If = quadva(@f3,[0,1],reltol,abstol);
fprintf('Approximate value of the integral = %f\n',If);
fprintf('True error        %e\n',abs(Itrue - If));
fprintf('%i (vector) evaluations of integrand were used.\n\n',nvfe);
disp('Using QUADL')
nvfe = 0;
If = quadl(@f3,0,1,abstol);
fprintf('Approximate value of the integral = %f\n',If);
fprintf('True error is      %e\n',abs(Itrue - If));
fprintf('%i (vector) evaluations of integrand were used.\n\n',nvfe);

%===Nested functions================================================
function y = f1(x)
  y = exp(-3*x) - cos(5*pi*x);
  nvfe = nvfe + 1;
end 

function y = f2(x)
  y = abs(x - 1/sqrt(3)) + abs(x + 1/sqrt(2));
  nvfe = nvfe + 1;
end

function y = f3(x)
  y = x.^(-2/3);
  nvfe = nvfe + 1;
end 
%===================================================================

end % Gladwellex