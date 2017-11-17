% Lista 3 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 3 - Ayiagari (1994)

clear all; clc; close all;
format short;

%% Calibracao
beta = 0.96;
alpha = 0.36;
u = @(c) (10*sqrt(c));
F = @(k,l) (k.^(alpha) .* l.^(1-alpha));
Fk = @(k,l) (alpha * (l./k).^(1-alpha)); % Pmg do capital
Fl = @(k,l) ((1-alpha) * (k./l).^(alpha)); % Pmg do trabalho
shocks = [0.1, 1];
n_shocks = length(shocks);
delta = 0.02;

pi_ll = 0.2; pi_lh = 1- pi_ll; % pi_ij: probabilidade de ocorrer j dado i
pi_hh = 0.6; pi_hl = 1 - pi_hh;

P = [pi_ll pi_lh;
     pi_hl pi_hh]; % Matriz de transição dos choques idiossincráticos
 
%% Trabalho do estado estacionario

% Encontrando distribuicao invariante dos choques idiossincraticos
epsilon = 1e-3;
dist_p = inv_dist_num(P, epsilon);

% Trabalho do steady state, considerando oferta de trabalho inelastica
L = dot(shocks, dist_p);
msg = strcat("Trabalho do steady state:  ", num2str(L));
disp(msg)

%% Espaco de estados
nk = 100;
k_grid = linspace(0, 10, nk);
k_linha_grid = linspace(0, 10, nk);

% Array de 3 dimensoes para variaveis de estado e a variavel de escolha
[k, neta, k_prime] = ndgrid(k_grid, shocks, k_linha_grid);

%% Resolvendo o modelo

K = 4;      % Chute inicial para capital do steady state
tol = 1e-4;
tol_v = 1e-3;   % tolerancia para value function
it = 0;
maxit = 2000;
maxit_v = 1000;     % maximo de iteracoes da value function
dist = tol + 1;
rho = 0.9 ;      % Parametro de estabilidade
V = zeros(nk, n_shocks);    % Chute inicial para funcao valor

tic
while dist > tol && it < maxit
    % Encontrando precos pelo problema da firma
    r = Fk(K, L);
    w = Fl(K,L);
    
    % Encontrando value function e policy function
    it_v = 0;   % contagem de iteracoes da value function
    dist_v = tol_v + 1;
    
    C = max((1 + r - delta) * k + w * neta - k_prime, 0);
    U = u(C);
    U([U == 0]) = NaN;     % Garantindo que o consumo nunca seja zero
    while dist_v > tol_v && it_v < maxit_v
        Y = permute(repmat(V, 1, 1, nk), [3, 2, 1]);
        EVh = beta * (Y(:, 1, :)*pi_hl + Y(:, 2, :)*pi_hh);
        EVl = beta * (Y(:, 1, :)*pi_ll + Y(:, 2, :)*pi_lh);
        H = U + [EVl, EVh];  % Possiveis candidatos a maximo
        [TV, ind] = max(H, [], 3);   % Procurando na terceira dimensao
        dist_v = norm(TV - V);
        it_v = it_v + 1;
        V = TV;
    end
    
    % Recuperando a funcao politica
    g = k_grid(ind);
    
    % Matriz de transicao no espaco produto
    
   gl = sparse(nk, nk);
   gh = sparse(nk, nk);
   
   for i = 1:nk
       gl(i, ind(i, 1)) = 1;
       gh(i, ind(i, 2)) = 1;
   end
   
   M = [gl * pi_ll, gl * pi_lh;
        gh * pi_hl, gh * pi_hh];
    
   M = full(M);
   
   % Encontrando a distribuicao invariante
   lambda = inv_dist_num(M, epsilon);
   
   % Capital do steady state
   Knew = dot([g(:,1)', g(:,2)'], lambda);
   
   % Atualizacao da iteracao
   it = it + 1;
   dist = abs(Knew - K);
   K = 0.1*Knew + 0.9*K;
   % K = Knew;
   
   if mod(it, 50) == 0
       disp(strcat("Iteracao: ", num2str(it)))
       disp(strcat("Capital: ", num2str(K)))
       disp(' ')
   end

end

disp('Fim do programa')
toc
