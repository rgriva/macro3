% Lista 1 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% Questao 1 - Modelo de McCall (1970)

clear all; close all; clc

%% Inicializando parametros e funcoes

beta = 0.98;
pi = 0.1;
b = 0;
wbar=20;
gamma=1/2;
alpha1 = 1/15;
alpha2 = -1/600;
u = @(c) (c.^gamma);
f = @(w) (alpha1 + alpha2.*w);

% Grid e parametros de iteracao
n = 1000;
wgrid = linspace(0,wbar,n);
h = wgrid(2) - wgrid(1);
epsilon = 1e-4;
maxit = 2000;
it = 0;
diff = 1;
V = zeros(n,1);
density = f(wgrid);

%% Iteracao - Item I
while diff > epsilon && it<maxit
    m = V'.*density;
    EV = (h/2)*(m(1) + m(n) + 2*sum(m(2:n-1)));
    N = u(b) + beta*EV;
    A = u(wgrid') + beta*((1-pi)*V + pi* V(1)*ones(size(V)));
    TV = max(N,A);
    ind = A > N;
    diff = sum(abs(TV - V));
    disp(diff);
    it = it + 1;
    V = TV;
end

g = ind;
wreserva = wgrid(max(find(g == 0)));

figure('DefaultAxesFontSize',16); 
subplot(2,1,1); 
plot(wgrid, V);
grid on;
ylabel('V(w)');
xlabel('w');
title('Funcao Valor')

subplot(2,1,2); 
plot(wgrid, g);
grid on;
ylabel('g(w)');
xlabel('w');
title('Funcao Politica')

%% Item II
alpha1 = 1/30;
alpha2 = 1/600;
f = @(w) (alpha1 + alpha2.*w);
density = f(wgrid);
it = 0;
diff = 1;
Vnew = zeros(n,1);

while diff > epsilon && it<maxit
    m = Vnew'.*density;
    EVnew = (h/2)*(m(1) + m(n) + 2*sum(m(2:n-1)));
    N = u(b) + beta*EVnew;
    A = u(wgrid') + beta*((1-pi)*Vnew + pi* Vnew(1)*ones(size(Vnew)));
    TVnew = max(N,A);
    ind = A > N;
    diff = sum(abs(TVnew - Vnew));
    disp(diff);
    it = it + 1;
    Vnew = TVnew;
end

gnew = ind;

figure('DefaultAxesFontSize',16); 
subplot(2,1,1); 
plot(wgrid, Vnew);
grid on;
ylabel('V(w)');
xlabel('w');
title('Funcao Valor - Nova Parametrizacao')

subplot(2,1,2); 
plot(wgrid, gnew);
grid on;
ylabel('g(w)');
xlabel('w');
title('Funcao Politica - Nova Parametrizacao')