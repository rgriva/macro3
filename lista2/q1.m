% Lista 2 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 1 - Bissecao e Newton - Raphson

clear all; clc; close all

%% Item A

f = @(x) (x.^(3) .* exp(-1 .* x.^(2)));
a = -10;
b = 9;
tol = 1e-10;
maxit = 1e7;

[root, it, value] = bissection_root(f, a, b, tol, maxit)