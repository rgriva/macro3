function [solution iterations]= gauss_jordan(A,b)
% GAUSS_JORDAN(A, b) uses the classic Gaussian method 
% for solving linear systems of the form Ax = b where A is nxn nonsingular 
% matrix and b is nx1.

E = [A, b];
[rows, cols] = size(A);
it = 0
% Triangulating the system
for j = 1:(cols-1)
    if E(j,j) ~= 0
        for i= (j+1):rows
            E(i,:) = E(i,:) - (E(i,j)/E(j,j))*E(j,:)
            it = it+1;
        end
    end
end

% Making backwards substitution
x = zeros(size(b));

x(length(b)) = E(length(b), cols+1)/E(rows,rows);
for k = length(b)-1:-1:1
    x(k) = (1/E(k,k)) * (E(k, cols+1) - E(k, k+1:cols)*x(k+1:length(b),1));
end

solution = x
iterations = it