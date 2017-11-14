% Lista 1 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 7 - Equity Premium Puzzle

clear all; close all; clc
tic;
%% Importando dados
% Dados estao na planilha data.xlsx, primeira sheet.

data = xlsread('data.xlsx',1);
consumo_bruto = data(:,2)*1e6; %A serie esta em milhoes de reais
ibov = data(:,3);
pop = data(:,4);
selic = data(:,5);
preco = data(:,7);

% Escolhendo de quando queremos comecar o exercicio numerico.
% Para que o ano base seja 1994, burn = 0.

burn = 0;

% Dividindo o consumo nominal em precos correntes pelo indice de precos
% e pela populacao
consumo_pc = consumo_bruto./(preco.*pop);
consumo_pc = consumo_pc(1+burn:end);

% Taxa de crescimento do consumo real per capta
gth_rate_c = (consumo_pc(2:end) - consumo_pc(1:end-1))./consumo_pc(1:end-1);

ibov_retorno = (ibov(2+burn:end) - ibov(1+burn:end-1))./ibov(1+burn:end-1);

%% Calibrando o modelo
% Seguindo a calibracao de Mehra & Prescott (1985)

mu = mean(gth_rate_c);
delta = std(gth_rate_c);
acf = autocorr(gth_rate_c,1);
phi = 0.5 * (1 + acf(2));
nbeta = 200;
nsigma = 200;
beta_grid = linspace(0.9,1,nbeta); 
beta_grid(nbeta) = 0.9999;
sigma_grid = linspace(1,10,nsigma);

x1 = 1 + mu + delta;
x2 = 1 + mu - delta;
x = [x1, x2];

% Ha 20 observacoes de retorno pois perdemos uma observacao ao calcular
% os retornos a partir dos indices. Por isso, vamos desconsiderar a SELIC
% de 1994, ano base.

diff_retorno = ibov_retorno - selic(2+burn:end);

% Premio de risco da economia
equity_premium = mean(diff_retorno);

%% Resolvendo o modelo
% Inicizalizando cadeia de Markov
phi11 = phi;
phi12 = 1- phi11;
phi21 = phi12;
phi22 = phi11;

P = [phi11, phi12; phi21 phi22];

% Encontrando a distribuicao invariante
[eigenvectors, eigenvalues] = eig(P, 'vector');
index = find(eigenvalues == 1);
inv_dist = eigenvectors(:,index)/sum(eigenvectors(:,index));

% Vetorizando as iteracoes
beta = repmat(beta_grid',1,length(sigma_grid));
sigma = repmat(sigma_grid, length(beta_grid), 1);
lambda1 = x1.^(-sigma);
lambda2 = x2.^(-sigma);

% Retornos livres de risco gerados pelo modelo, para cada combinacao
% de beta (linha) e sigma (coluna).

p1f = (P(1,1)*lambda1 + P(1,2)*lambda2).*beta;
R1f = 1./p1f - 1;
p2f = (P(2,1)*lambda1 + P(2,2)*lambda2).*beta;
R2f = 1./p2f - 1;

% Usando a distribuicao invariante da cadeia de Markov
Rf = inv_dist(1)*R1f + inv_dist(2)*R2f;

% Calculando os pesos dos estados
Re = zeros(size(Rf));
A = zeros(length(x),length(x));
b = zeros(length(x),1);
R = zeros(size(A));

% Encontrando valores de Re para cada combinacao de beta e sigma
for k = 1:nbeta
    for m = 1:nsigma
        for i = 1:length(x)
            for j = 1:length(x)
                A(i,j) = beta_grid(k)*(x(j)^(1-sigma_grid(m))*P(i,j));
            end
            b(i) = sum(A(i,:));
        end
        q = (eye(size(A)) - A)\b;
        Re11 = (q(1) + 1)*x(1)/q(1) - 1;
        Re12 = (q(2) + 1)*x(2)/q(1) - 1;
        Re21 = (q(1) + 1)*x(1)/q(2) - 1;
        Re22 = (q(2) + 1)*x(2)/q(2) - 1;
        R = [Re11 Re12; 
             Re21 Re22];
        Re1 = P(1,:)*R(1,:)';
        Re2 = P(2,:)*R(2,:)';
     
        Re(k,m) = inv_dist(1)*Re1 + inv_dist(2)*Re2;
        
    end
end
toc

%% Analisando a simulacao
equity_premium_simulado = Re - Rf;

dist = abs((equity_premium_simulado - equity_premium*ones(size(Rf))));
[a,b] = find(dist == min(dist(:)));
beta_optimal = beta_grid(b)
sigma_optimal = sigma_grid(a)

%% Plotando
figure('DefaultAxesFontSize',16); hold on
plot = surf(beta_grid,sigma_grid,dist);
grid on;
set(plot,'LineStyle','none')
title('Equity Premium')
xlabel('Beta')
ylabel('Sigma')
zlabel('Modulo da Diferenca')