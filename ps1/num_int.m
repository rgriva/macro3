function integral = num_int(f, a, b, steps)
% NUM_INT(f, a, b, steps) computes the numerical integral of a real 
% function f between a and b and the
% number of steps for the discretization of the domain, 
% using the Trapezoidal rule.

grid = linspace(a,b,steps);
h = grid(2) - grid(1);
v = f(grid);
integral = (h/2)*(sum(v) + sum(v(2:(steps -1))));