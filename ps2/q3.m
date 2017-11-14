% Lista 2 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 3 - Wallace JME 1992

clear all; clc; close all

%% Item A

% Parametros
u = @(x) sqrt(x);
v = @(x) sqrt(x);
w=1;
finv = @(y) (-2.*y.^2 + 2.*y.*sqrt(w + y.^2))
beta = 0;
delta = 0.05;
I = 2;
J = input('Defina o numero de estados de xt (J): ');
lambda = input('Defina o parametro lambda: ');

pi = 1/I;
theta = 1/J;

% Calculando os estados
gridI = 1:I;
gridJ = 1:J;
lnN = beta + delta*gridI;
N = exp(lnN);
lnx = lambda + delta*gridJ;

K = ((1-I):(J-1))';

Z = exp(lambda - beta + delta*K); % vetor coluna

%% Calculando probabilidades
probx = theta*ones(J,1);
probN = theta*ones(I,1);
phi_k = zeros(size(K));     % P(Zt = Zk)



p = zeros(J,1);
for m = 1:length(K)
    for j=1:J
        if j-K(m) > 0 && j-K(m) < I + 1
            p(j,1) = theta * pi;
        else
            p(j,1) = 0;
        end
    end
    phi_k(m) = sum(p);
end

% P(Nt = Ni | Zt = Zk), Z nas linhas e N nas colunas
phi_ki = zeros(length(K),I);

for i = 1:I
    for m = 1:length(K)
        if i+K(m) > 0.5 && i+K(m) < J+0.5 
        phi_ki(m,i)=(pi*theta)/phi_k(m);
        end
    end
end

%% Iteracao da funcao valor

tol = 1e-5;
it = 0;
dist = 10;
maxit = 1e4;
y0 = ones(length(K),1);
phizao = repmat(phi_k', 2,1);

while dist > tol && it < maxit
    M = y0./Z;
    m = repmat(M',I,1);
    n = repmat(N',1,length(M));
    g = 0.5*sqrt(m./n); % cada linha é avaliada pra um N e cada coluna pra um h
    b = sum(phizao.*g, 2); % somando nas linhas
    Eg = phi_ki * b;
    Ty0 = finv(Eg);
    
    dist = norm(Ty0 - y0);
    it = it + 1;
    y0 = Ty0;
end

y0;

%% Plotando z/y(z)

figure('DefaultAxesFontSize',16); hold on; grid on;
scatter(K, 1./M, [] ,'r', 'filled')
xlabel('z')
ylabel('z/y(z)')
title('J = 40, lambda = 0')