% Lista 1 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 3 - Crescimento Estocastico

clear all; close all; clc

%% Inicializando parametros e funcoes
alpha = 0.7;
beta = 0.98;
gamma = 2;
delta = 0.1;
zh = 1.2;
zl = 0.8;
xi = 0.7; % Prob(z = zh | zh)
zeta = 0.2; % Prob(z = zl | zl)
f = @(k,z) (z.*k.^alpha);
u = @(c) ((c.^(1-gamma))/(1-gamma));

% Grids
kssh = (1/(zh*alpha)*(1/beta-1+delta))^(1/(alpha-1));
kssl = (1/(zl*alpha)*(1/beta-1+delta))^(1/(alpha-1));
n = 1000;
kgrid = linspace(0.3*kssl, 2*kssh, n);
kprimegrid = kgrid;

% Parametros de iteracao
epsilon = 1e-4;
maxit = 2000;
it = 0;
diff = 1;
Vh = zeros(n,1); % Primeira coluna se refere a Low e a segunda a High
Vl = zeros(n,1);
%% Iteracao
tic;
k = repmat(kgrid', 1, n);
kprime = repmat(kprimegrid, n, 1);
Ch = max((1-delta)*k + f(k,zh) - kprime, 0);
Cl = max((1-delta)*k + f(k,zl) - kprime, 0);
Uh = u(Ch);
Ul = u(Cl);

while diff > epsilon && it < maxit
    Rl = Ul + beta*(zeta*repmat(Vl', n,1) + (1-zeta)*repmat(Vh',n,1));
    Rh = Uh + beta*((1-xi)*repmat(Vl', n,1) + xi*repmat(Vh',n,1));
    [TVl, indl] = max(Rl,[],2); 
    [TVh, indh] = max(Rh,[],2);
    diff = sum(abs(TVl - Vl)) + sum(abs(TVh - Vh));
    disp(diff)
    it = it +1;
    Vl = TVl;
    Vh = TVh;
end
toc 

%Recuperando funcao politica
gh = k(indh);
gl = k(indl);

%% Plotando resultados
figure('DefaultAxesFontSize',16);
subplot(2,1,1); 
plot(kgrid, Vh);
hold on
plot(kgrid, Vl);
ylabel('V(k,z)');
xlabel('k');
title('Funcao Valor')
legend('V(k, z alto)', 'V(k, z baixo)', 'Location', 'SouthEast')

subplot(2,1,2); 
plot(kgrid, gh);
hold on;
plot(kgrid, gl);
ylabel('g(k,z)');
xlabel('k');
title('Funcao Politica')
legend('g(k, z alto)', 'g(k, z baixo)', 'Location', 'SouthEast')