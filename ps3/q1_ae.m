% Lista 3 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 1 - Hugget (JEDC - 1993)
% Itens (a)-(e)

clear all; clc; close all

%% Calibracao
eh = 1;
el = 0.1;
shocks = [el, eh];
pi_hh = 0.925;      % pi_ji: probabilidade de ocorrer i dado j.
pi_hl = 1 - pi_hh;
pi_lh = 0.5;
pi_ll = 1 - pi_lh;
P = [pi_hh, pi_hl;
     pi_lh, pi_ll;];
n_shocks = length(shocks);
% clearvars pi_hh pi_hl pi_lh pi_ll
 
beta = 0.99322;
sigma = 1.5;
%a_bar = input('Escolha um valor de a_barra (limite de credito): ');
%q = input('Escolha um valor inicial de q (positivo): ');  % Chute inicial do preco q
a_bar = -2;
q = 1.0124;
u = @(c) (c.^(1 - sigma))/(1- sigma);
n = 100;
a_grid = linspace(a_bar, -3*a_bar, n);

%% Iterando funcao valor
[a, e, a_linha] = ndgrid(a_grid, shocks, a_grid); % Discretizando o espaco de estados 
epsilon = 1e-3;
maxit = 2000;
it = 0;
dist = epsilon + 1;

C = max(a + e - a_linha*q,0); 
U = u(C);
V = zeros(n, n_shocks);

tic
while dist > epsilon && it < maxit
    Y = repmat(V, 1, 1, n);
    J = permute(Y, [3, 2, 1]); % Colocando na dimensao correta
    EVl = beta*(J(:,1,:)*pi_ll + J(:,2,:)*pi_lh);
    EVh = beta*(J(:,1,:)*pi_hl + J(:,2,:)*pi_hh);
    EV = [EVl, EVh];
    H = U + EV; % Candidatos a maximo
    [TV , index] = max(H, [], 3);
    
    dist = norm(TV - V);
    it = it + 1;
    V = TV;
end
g = a_grid(index);
toc

%% Plotando
figure('DefaultAxesFontSize',16);
subplot(2,1,1);
plot(a_grid, V(:,1))
hold on; grid on;
plot(a_grid,V(:,2))
title('Funcao valor')
legend('e baixo', 'e alto', 'Location', 'SouthEast')
xlabel('Ativo a')

subplot(2,1,2); grid on
plot(a_grid, g(:,1))
hold on; grid on;
plot(a_grid, g(:,2))
title('Funcao politica')
legend('e baixo', 'e alto', 'Location', 'SouthEast')
xlabel('Ativo a')

%% Matriz de transicao M

gl = sparse(n, n);
gh = sparse(n, n);


for i=1:n % Gerando funcoes indicadoras
    gl(i, index(i,1)) = 1;
    gh(i, index(i,2)) = 1;
end
M = zeros(2*n);

M = [gl*pi_ll gl*pi_lh; gh*pi_hl gh*pi_hh];
M = full(M);

%% Autovalores e autovetores de M

% Metodo do maior autovalor
[eig_vec, eig_values] = eig(M');
[r,c] = find(eig_values == max(max(eig_values)));
max_eigenvalue = eig_values(r,c)
inv_dist = eig_vec(:,r)/sum(eig_vec(:,r))

% Metodo da iteracao
tol = 1e-10;
init = rand(length(M),1);
x0 = init/sum(init);
dif = tol + 1;
it = 0;

while dif > tol && it < maxit
    new = M'*x0;
    dif = max(abs(new - x0));
    it = it + 1;
    x0 = new;
end

inv_dist_it = new/sum(new)

%% Excesso de credito

B = [g(:,1)' g(:,2)']*inv_dist_it    
    
    
    
    
    
    
    
    
    
    
    
