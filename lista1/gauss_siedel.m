function [solution, iterations] = gauss_siedel(A, b, x0, epsilon, maxit)
% GAUSS_SIEDEL(A, b, x0, epsilon, maxit) uses the classic Gauss-Siedel method 
% for solving linear systems of the form Ax = b where A is nxn nonsingular 
% matrix, b is nx1, x0 is the initial guess for the solution, 
% epsilon is a convergence parameter and maxit is the maximal 
% number of iterations.
% [x, y] = GAUSS_SOLVER(A, b, x0, epsilon, maxit) returns the computed 
% solution at x and the total number of iterations at y.
% Sufficient conditions for convergence: Matrix A should be diagonal 
% dominant or positive semi-definite.

L = tril(A);
U = A - L;

diff = 10;
it = 0;
xold = x0;

while diff > epsilon && it < maxit
    xnew = inv(L)*(b - U*xold)
    diff = sum(abs(xnew - xold));
    it = it + 1;
    xold = xnew;
end

solution = xnew;
iterations = it;