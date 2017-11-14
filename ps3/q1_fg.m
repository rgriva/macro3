% Lista 3 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 1 - Hugget (JEDC - 1993)
% Itens (f)-(g)

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

beta = 0.99322;
sigma = 1.5;

a_bar = input('Limite inferior do credito: ');
qh = 1.02; % Chutes iniciais para limites superior e inferior de q
ql = 0.99;
q = ql;
u = @(c) (c.^(1 - sigma))/(1- sigma);
n = 100;
a_grid = linspace(a_bar, -3*a_bar, n);

[a, e, a_linha] = ndgrid(a_grid, shocks, a_grid); % Discretizando o espaco de estados 
epsilon = 1e-5;

B = 1;
m = 0.0025;
error = 1 + m;
max_it = 2000;
iterada = 0;
V = zeros(n, n_shocks);

%% Encontrando preco de equilibrio

while abs(B) > m && iterada < max_it
    C = max(a + e - a_linha*q,1e-18); 
    U = u(C);
    dist = epsilon + 1;
    it = 0;
        while dist > epsilon && it < max_it
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
        
    gl = sparse(n, n);
    gh = sparse(n, n);


    for i=1:n % Gerando funcoes indicadoras
        gl(i, index(i,1)) = 1;
        gh(i, index(i,2)) = 1;
    end

    M = [gl*pi_ll gl*pi_lh; gh*pi_hl gh*pi_hh];
    M = full(M);
    
    [eig_vec, eig_values] = eig(M');
    [r,c] = find(eig_values == max(max(eig_values)));
    max_eigenvalue = eig_values(r,c);
    inv_dist = eig_vec(:,r)/sum(eig_vec(:,r));
    
    B = [g(:,1)' g(:,2)']*inv_dist;

    if B > 0
        ql = ql + 0.5*(q - ql);
    else
        qh = qh - 0.5*(qh - q);
    end
    
    q = 0.5*(qh + ql);
    iterada = iterada + 1;
    
    if mod(iterada,10) == 0
        disp(iterada)
        disp(q)
        disp(B)
    end
end
    
    
    
    
    
    
    
    
    
    
    