function [solution, iterations] = jacobi_solver(A, b, x0, epsilon, maxit)
% JACOBI_SOLVER(A, b, x0, epsilon, maxit) uses the classic Jacobi's method 
% for solving linear systems of the form Ax = b where A is nxn nonsingular 
% matrix, b is nx1, x0 is the initial guess for the solution, 
% epsilon is a convergence parameter and maxit is the maximal 
% number of iterations.
% [x, y] = JACOBI_SOLVER(A, b, x0, epsilon, maxit) returns the computed 
% solution at x and the total number of iterations at y.
% ATTENTION: Matrix A should be diagonal dominant in order to this work 
% properly regardless of the initial condition.

D = diag(diag(A));
R = A - D;

diff = 10;
it = 0;
xold = x0;

while diff > epsilon && it < maxit
    xnew = inv(D)*(b - R*xold)
    diff = sum(abs(xnew - xold));
    it = it + 1;
    xold = xnew;
end

solution = xnew
iterations = it