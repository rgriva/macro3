% Lista 1 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 5 - Sistemas Lineares

clear all; close all; clc

%% Item A

A = [5, -2, 3; -3, 9, 1; 2 -1 -7];
b = [-1;2;3];
x0 = zeros(3,1);
maxit = 1000;
epsilon = 1e-3;

[sol_jacobi, it_jacobi] = jacobi_solver(A,b,x0,epsilon,maxit)

%% Item B
[sol_gauss, it_gauss] = gauss_jordan(A,b)