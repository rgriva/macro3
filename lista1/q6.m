% Lista 1 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 6 - Sistemas Lineares

clear all; close all; clc

%% Item A

A = [1, -5; 7, -1];
b = [-4;6];
x0 = zeros(2,1);
maxit = 1000;
epsilon = 1e-3;

[sol_jacobi, it_jacobi] = jacobi_solver(A,b,x0,epsilon,maxit)

%% Item B
x0 = [-3 ; 1];
[sol_jacobi, it_jacobi] = jacobi_solver(A,b,x0,epsilon,maxit)

%% Item C
[sol_gauss, it_gauss] = gauss_jordan(A,b)