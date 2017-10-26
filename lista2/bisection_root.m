function [root, it, value] = bisection_root(f, a, b, tol, maxit)
% BISSECTION_ROOT uses bissection method to finds the root of a 
% univariate continuous @-function f in [a, b], provided that exists one, 
% with tol being the tolarance parameter for convergence and maxit a 
% maximal number of iterations.

% Checking if the method works.
if sign(f(a)) == sign(f(b))
    message = 'Root not assured to exist';
    error(message)
    root = NaN;
end

% Checking extremum points
if abs(f(a)) < tol
    root = a;
    it = 0;
elseif abs(f(b)) < tol
    root = b;
    it = 0;
end

% Starting iteration
it = 0;
dist = 1;

while dist > tol & it < maxit
    xm = 0.5 * (a + b);
    if sign(f(xm))*f(b) < 0
        a = xm;
    else
        b = xm;
    end
    dist = abs(f(xm));
    it = it + 1;
end

% root and final value attained
root = xm;
value = f(xm);




