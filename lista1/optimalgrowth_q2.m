% Lista 1 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 2 - Crescimento Neoclassico

clear all; close all; clc

%% Inicializando parametros e funcoes
alpha = 0.7;
beta = 0.98;
gamma = 2;
delta = 0.1;
f = @(k) (k.^alpha);
u = @(c) ((c.^(1-gamma))/(1-gamma));
kss = (1/(alpha)*(1/beta-1+delta))^(1/(alpha-1)); % Resultado analitico
 
% Criando grids
n = 5000;
kgrid = linspace(0.01*kss, 1.25*kss, n);
kprimegrid = kgrid;

% Parametros de iteracao
epsilon = 1e-4;
maxit = 2000;
it = 0;
diff = 1;
V = zeros(n,1); % Chute inicial para V e para sua imagem por T
TV = zeros(n,1);
index = zeros(n,1); % Indexador de kgrid para g

%% Iterando a funcao valor
tic;
% Calculando consumos
K = repmat(kgrid',1,n);
Kprime = repmat(kgrid, n,1);

C = max(f(K) + (1-delta)*K - Kprime, 0); % Todos os consumos possiveis
U= u(C); % Todas as utilidades possiveis

while diff > epsilon && it < maxit
    D = U + beta*repmat(V', n, 1);
    [TV, index] = max(D, [], 2); % Procurando o maximo linha por linha
    diff = max(abs(TV - V));
    disp(diff)
    it = it + 1;
    V = TV;
end

g = kgrid(index);
toc
%% Plotando resultados
figure; 
subplot(2,1,1); 
plot(kgrid, V);
grid on;
ylabel('V(k)');
xlabel('k');
title('Funcao Valor')

subplot(2,1,2); 
plot(kgrid, g);
grid on;
ylabel('g(k)');
xlabel('k');
title('Funcao Politica')

%% Simulando trajetoria do capital
k0 = 2;
d = abs(kgrid - k0*ones(size(kgrid)));
start = find(d == min(d));

knew = g(start);
it = 1;
path = zeros(maxit, 1);
path(1) = knew;
while abs(knew - k0) > 0.001 && it < maxit
    k0 = knew;
    knew = g(find(kgrid == k0));
    it = it + 1;
    path(it) = knew;
end

path = path(1:it);
figure; plot(path)
grid on;
title('Caminho do Capital - Simulado')
xlabel('Tempo')
ylabel('Estoque de Capital')

kss_numerical = path(it)
css_numerical = f(kss_numerical) + (1-delta)*kss_numerical - kss_numerical
yss_numerical = f(kss_numerical)

%% Dependencia de beta
betagrid = linspace(0, 1, 100);
kss_beta = (1/(alpha)*(1./betagrid-1+delta)).^(1/(alpha-1));
figure; plot(betagrid, kss_beta)
grid on
xlabel('Beta')
ylabel('Capital do Estado Estacionario')
title('Relacao Kss e Beta')

%% Dependencia de gamma
gammagrid = linspace(0, 10, 100);
kss_beta = (1/(alpha)*(1./beta-1+delta)).^(1/(alpha-1))*ones(size(gammagrid));
figure; plot(gammagrid, kss_beta)
grid on
xlabel('Gamma')
ylabel('Capital do Estado Estacionario')
title('Relacao Kss e Gamma')