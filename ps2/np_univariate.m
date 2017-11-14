function [root, it, value] = np_univariate(f, a, b, tol, maxit)
% NP_UNIVARIATE uses Newton-Raphson method to try to find a root of a 
% univariate continuous symbolic function f in [a, b], provided that there 
% exists one,  with tol being the tolarance parameter for convergence and 
% maxit a  maximal number of iterations.

% Random initial condition
x0 = 0.5*(a+b);

flinha = diff(f);
f_handle = matlabFunction(f);
flinha_handle = matlabFunction(flinha);

if f_handle(a)*f_handle(b) > 0
    message = 'Root not assured to exist';
    error(message)
    root = NaN;
    it = 0
    value = NaN;
end

it = 0;
dist = 1;

while dist > tol & it < maxit
    x1 = x0 - (f_handle(x0))/(flinha_handle(x0));
    
    % Checking if test is not out of bounds and reshuffling it
    if x1 < a || x1 > b
        x1 = a + (b-a)*rand;
    end
    
    deltax = abs(x1-x0);
    deltaf = abs(f_handle(x1));
    dist = max(deltax, deltaf);
    it = it + 1;
    x0 = x1;
end

root = x1;
value = f_handle(root);
    