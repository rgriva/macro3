% Lista 2 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 2 - Newton-Raphson multivariado

% ATENCAO: Fazemos forte uso da toolbox de Matemática Simbólica do Matlab!

clear all; clc; close all

%% Definindo funcoes

syms x y
f(x,y) = 4/((x-1)^2 + 4*y^2 + 1)
f_handle = matlabFunction(f);
Dfx (x,y) = diff(f,x);
Dfy (x,y) = diff(f,y);
jacobian (x,y) = [Dfx(x,y); Dfy(x,y)];

Dfxx = diff(Dfx, x);
Dfxy = diff(Dfx, y);
Dfyy = diff(Dfy, y);
hessian (x,y) = [Dfxx(x,y), Dfxy(x,y);
                 Dfxy(x,y), Dfyy(x,y)];
hinv (x,y) = inv(hessian)

%% Operadores numericos
% Lidar com operadores numericos ao inves de simbolicos aumenta a
% performance do codigo

Hinv = matlabFunction(hinv);
gradiente = matlabFunction(jacobian);

%% Plotando f para analisar visualmente
gridx = linspace(-10,10,100);
gridy = linspace(-10,10,100);
[x,y] = meshgrid(gridx, gridy);
figure; hold on; grid on
title('F(x)')
xlabel('x')
ylabel('y')
data = f_handle(x,y);
surf(gridx,gridy,data)
             

%% Iteracao

it = 1;
dist = 10;
tol = 1e-7;
maxit = 1e5;

% Chute inicial, vetor coluna
xold = [0, 0]';
xk = zeros(maxit,1); xk(1) = xold(1);
yk = zeros(maxit,1); yk(1) = xold(2);
fk = zeros(maxit,1); fk(1) = f(xold(1),xold(2));

% Definindo bounds
u_boundx = 5;
l_boundx = -5;
u_boundy = 5;
l_boundy = -5;

tic
while dist > tol && it < maxit
    xnew = xold - Hinv(xold(1), xold(2)) * gradiente(xold(1), xold(2));
    
    % Garantindo que estamos dentro de um compacto em R2:
    if xnew(1) > u_boundx || xnew(1)< l_boundx
        xnew(1) = l_boundx + (u_boundx - l_boundx)*rand;
    elseif xnew(2) > u_boundy || xnew(2) < l_boundy
        xnew(2) = l_boundy + (u_boundy - l_boundy)*rand;
    end
        
    it = it + 1;
    xk(it) = xnew(1);
    yk(it) = xnew(2);
    fk(it) = f_handle(xk(it), yk(it));
    dist = norm(xnew - xold);
    xold = xnew;
end

toc

xk = xk(1:it);
yk = yk(1:it);
fk = fk(1:it);
argmax = xnew
value = fk(it)
it
    
