function [root, it, value] = bissection_root(f, a, b, tol, maxit)
% BISSECTION_ROOT uses bissection method to finds the root of a 
% univariate continuous @-function f in [a, b], providaded that exists one, 
% with tol being the tolarance parameter for convergence and maxit a 
% maximal number of iterations.

% Checking if the method works.
if sign(f(a)) == sign(f(b))
    message = 'Root not assured to exist';
    disp(message)
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
c = 0.5 * (a + b);
while abs(f(c)) > tol & it < maxit
    if sign(f(c)) ~= sign(f(b))
        c = 0.5 * (c + b)
    else
        c = 0.5 * (c + a)
    end
    it = it + 1
end

% root and final value attained
root = c;
value = f(c);




