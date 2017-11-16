% Lista 3 - Macroeconomia III 2017
% Alunos: Alexandre Machado e Raul Guarini
% (Fortemente baseado no trabalho de Katia Alves e Ricardo Elias)
% Questao 2 - Imrohoroglu (1992)

clear all; clc; close all

%% Parametros

theta  = 0.25;      
sigma = 1.5;      
beta   = 0.995;
pi=input('Escolha pi =  ');     

P   = [0.9565,0.0435;0.5,0.5];     

tau = pi;                 % tx. de crescimento da base monetária

inc = 0.027;               % incremento 
g   = 301;                 % tamanho do grid   

tol = 1e-6;                         
R   = 1/(1+pi);                    

% fixando a oferta de moeda
M       = 5;      

% Discretização do espaço das v.estado
m = linspace(0,inc*300,g);
m = m';
m_hj = m*ones(1,g);
m_novo = m_hj';

y = ones(g,g);
one = ones(g,1);
dist = 100;
it = 1;

%% Resolvendo o modelo
while (it <= 100) & (dist > tol)
   
    MM(it)       = M;
    M_med      = M*ones(g,g);
   
    cemp          = R*m_hj + y + tau*R*M_med - m_novo;
    u_emp      = (cemp.^(1-sigma) - 1) / (1-sigma);
    j   = [];   k = [];
    [j,k]       = find(cemp <= 0);
    siz         = size(j);
   
    for r = 1:siz(1,1)
        u_emp(j(r,1),k(r,1)) = -inf;
    end
   
    cdes          = R*m_hj + theta*y + tau*R*M_med - m_novo;
    u_desemp     = (cdes.^(1-sigma) - 1) / (1-sigma);   
    j   = [];   k = [];
    [j,k]       = find(cdes <= 0);
    siz         = size(j);
   
    for r = 1:siz(1,1)
        u_desemp(j(r,1),k(r,1)) = -inf;
    end
       
    v_emp  = zeros(1,g);
    v_desemp = zeros(1,g);
    pos     = zeros(g,2);
    dist1   = 100;
    dist2   = 100;
 
    while (dist1 > .01) | (dist2 ~= 0)
       
        vala = (u_emp + beta*one*( P(1,1)*v_emp + P(1,2)*v_desemp ) )';
        valb = (u_desemp + beta*one*( P(2,1)*v_emp + P(2,2)*v_desemp ) )';
       
        [V_empl,pos1]   = max( vala );
        [V_unemp,pos2]  = max( valb );
       
        V           = [V_empl,V_unemp];
        v           = [v_emp,v_desemp];
        Pos         = [pos1',pos2'];
       
        dist1       = abs(V - v);
        dist1       = max( (max(dist1))' );
        dist2       = max(any(Pos-pos));
       
        v_emp      = V_empl;
        v_desemp     = V_unemp;
        pos         = Pos;
       
    end
       
    posa    = zeros(g,g);
    posb    = zeros(g,g);
 
    for i = 1:g
        posa(i,pos1(1,i)) = 1;
        posb(i,pos2(1,i)) = 1;
    end
   
    tran_mat    = [P(1,1)*posa,P(1,2)*posa;P(2,1)*posb,P(2,2)*posb];
    lambda           = (ones(2*g,1))/(2*g);
 
    dist1   = 100;                  %achando a matriz de dist invariante
    while dist1 > tol
        lambda_p     = tran_mat'*lambda;
        dist1   = max(abs(lambda_p-lambda));
        lambda       = lambda_p;
    end
   
    M_p     = [m;m]'*lambda;
    dist    = abs((M_p-M)/M);
    M       = M_p;
    it       = 1 + it;
   
end
 
 
one     = ones(2*g,1);
m_med  = M;
m_var   = lambda' * (([m;m]-M*one).^2);
m_dp   = sqrt(m_var);

%% Estatisticas do Modelo
for i=1:g
    % consumption
    consemp(i,1)   =cemp(i,pos(i,1));
    consdes(i,1)   =cdes(i,pos(i,2));
    % utility
    util_emp(i,1) = u_emp(i,pos(i,1));
    util_des(i,1) = u_desemp(i,pos(i,2));
end

cons=[consemp;consdes];
c_med  = lambda'*cons;
c_dp   = sqrt(lambda' * ((cons-c_med*one).^2));

yy      = [ones(g,1);.25*ones(g,1)];
y_med  = lambda' * yy;
y_dp   = sqrt(lambda' * ((yy-y_med*one).^2));

u   = [util_emp;util_des];
u_med  = lambda' * u;

%% Tabela
disp(' ')
info = sprintf('    real cash balances médio          =     %3.4f', m_med);
disp(info);
disp(' ')
info = sprintf('    dp do real cash balance           =     %3.4f', m_dp);
disp(info);
disp(' ')
info = sprintf('    consumo médio                     =     %3.4f', c_med);
disp(info);
disp(' ')
info = sprintf('    dp do consumo                     =     %3.4f', c_dp);
disp(info);
disp(' ')
info = sprintf('    renda média                       =     %3.4f', y_med);
disp(info);
disp(' ')
info = sprintf('    dp da renda                       =     %3.4f', y_dp);
disp(info);
disp(' ')
info = sprintf('    utilidade média                   =     %3.4f', u_med);
disp(info); disp(' ')
 
 
inc         = m(2,1)-m(1,1);
pol_emp    = inc*( pos(:,1) )';
pol_desemp   = inc*( pos(:,2) )';
mp          = m';
 
Dist        = ( lambda(1:g,1) )' + ( lambda(g+1:2*g,1) )';
