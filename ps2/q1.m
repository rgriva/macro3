% Lista 2 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 1 - Bissecao e Newton - Raphson

clear all; clc; close all

%% Item A

syms x
f = x^(3) * exp(-x^(2));
flinha = matlabFunction(diff(f));


% Procurando no lugar certo
a = 0.5;
b = 6;
tol = 1e-6;
maxit = 1e7;

tic
[root, it, value] = bisection_root(flinha, a, b, tol, maxit);

f_handle = matlabFunction(f);

max_global = f_handle(root)
argmax = root
it
toc

%% Item B

syms x
f = x^(3) * exp(-x^(2));
flinha = diff(f);

a = 0.5;
b = 6;

tic
[root, it, value] = np_univariate(flinha, a, b, tol, maxit)

f_handle = matlabFunction(f);

max_global = f_handle(root)
argmax = root
it
toc
