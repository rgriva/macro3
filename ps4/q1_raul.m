% Lista 4 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 1 - Diamond and Dybvig (Matriz F)

clc; clear all; close all

%% Calibração

N = 2;
R = 1.05;
delta = 2;
u = @(c) (c.^(delta))/(1-delta);
e = 1; % dotacao per capta
p = 0.5;

j_set = 0:1:N;  % conjunto dos possiveis valores de |w| 

%% Item A
F = zeros(N+1, N+1);    

F(:, N+1) = j_set*R^((1-delta)/delta);  % coeficientes da ultima data
for m = 1:N
    for k = m:N
        F(m, k) = (p * (1+F(m, k+1))^(delta) + (1-p) * F(m+1, k+1)^(delta))^(1/delta);
    end
end
disp('Matriz F')
disp(F)

%% Item B

% Calculando pagamento otimo na primeira data:
% Linha: quantidade de pacientes ate aquela posicao na fla
% Coluna: posicao na fila

% Sorteando um estado da natureza
omega = binornd(1, p, N, 1)';
disp('Fila: ')
disp(omega)

% Coletando os anuncios
cont_pac = zeros(1, N+1);
for i = 1:N
    cont_pac(i) = sum(omega(1:i));
end
cont_pac(end) = cont_pac(end - 1); % ultimo termo representa a ultima data
cont_pac

z = N*e;
x = zeros(1,N);

for i=1:N
    x(i) = z/(1 + F(1 + cont_pac(i-1) , i));
    z = z - x(i);
end

disp('Pagamento otimo para primeira data')
x

% Calculando pagamento otimo na ultima data, igualitario
y = R*z/cont_pac(end)