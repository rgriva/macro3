clear all;
close all;
clc;

%% QUESTAO 1

% Dados
A=1;
N=2;
R=1.05;
delta=2;
P=0.5;
a=1.11;
e=1;

% Utilidade
u=@(c)((c.^(1-delta))/(1-delta));

% Chute inicial para o multiplicador 

Lambda=10;

dist=1;
tol=10^(-6);

itmax=10^(6);
it=0;

while it<itmax && dist>tol
    
      it=it+1;
      
      gamma=abs(((1+Lambda/(1-P))/(A-Lambda/P))^(1/delta));
      
      M=N+1;
      F=zeros(M,M);
      
    % Obs:. Para evitar confusão com os índices, seguiremos o padrão
    % estabelecido pelo handout, utilizando "i" para representar
    % o agente, e "j" para representar o seu anuncio.
    
    % Aqui vamos calcular o F ótimo para o problema do agente após a
    % realização de um certo (w)
    
    for j=1:M
        F(j,M)=gamma*(j-1)*R^(1/delta-1);
    end
    
    % Agora vamos completar os temos da matriz F, a ideia e resolver o
    % problema retroativamente utilizando a informação calculada na parte
    % anterior.
      
      for k=1:M-1;
        j=M-k;
        for i=1:M-k
            F(i,j)=(P*((F(i,j+1)+1)^delta)+(1-P)*(F(i+1,j+1)^delta))^(1/delta);
        end    
      end
      
          
    % Precisamos ver qual seria o ganho do agente de desviar para, então,
    % podermos atualizar o multiplicador e iterar até chegarmos ao ótimo:
      
      Ganho=zeros(M,M); 

    % Vamos fazer da mesma forma como fizemos para calcular F. Aqui o ganho
    % de desviar para o último indivíduo é equivalente ao que ele ganharia
    % em t=1 se desviasse do mecanismo calculado acima.
    
    for j=1:M 
        Ganho(j,M) = (P/(1-P))*(j-1)^delta*R^(1-delta);
    end
      
    % Completando a matriz:  
    
      for k=1:M-1 
          j=M-k;
        for i=1:M-k
            Aux=Ganho(i,j+1)*(F(i,j+1)^(1-delta));
            if isnan(Aux)==1
                Aux=0;
            end    
            Ganho(i,j)=P*(Aux-1)*(1+F(i,j+1)^(delta-1))+(1-P)*(Ganho(i+1,j+1));
        end    
      end
      
    % Agora é só calcular o novo lambda, seguindo a metodologia da função
    % penalidade, e verificar a sua convergência (esse loop maior)
    
    lambda=exp(Ganho(1,1))*Lambda;
    dist=abs(Lambda-lambda);
    Lambda=lambda; 
end


%% Mecanismo otimo

% A matriz abaixo representa todo o espaço de "estados" possíveis.

Omega=ones(length(N));

for k=1:N
    Omega(k)=round(rand(1));
end 

NewF=zeros(N,N); 
NewF=F(1:N,2:N+1);
 
% Tx. Poupança
S=ones(1,N);
 
% Vamos verificar agora quantos agentes pacientes temos na economia:

Pacientes(1:N)=zeros(1,N); 

for i=1:N-1
    Pacientes(i+1)=sum(Omega(1:i)); 
end

% Vamos verificar agora, dada nossa matriz Omega gerada anteriormente,
% onde estão os agentes impacientes;

Impacientes=find(Omega==0);

% Poupança

for i=Impacientes
    S(i)=NewF(Pacientes(i)+1, i);
    S(i)=S(i)/(1+S(i));
end

% Logo, S representa a taxa de poupança do banco.

Acumulo=1;
Dotacao_Economia=e*N;
Consumo=Dotacao_Economia*ones(1,N);

for i=1:N
    Acumulo=Acumulo*S(i);
    Consumo(i)=Acumulo*Dotacao_Economia;
end

Pag_Impacientes=(1-S).*[Dotacao_Economia,Consumo(1:N-1)];

Pag_Pacientes=zeros(1,N);
if sum(Omega)>0
    Aux=Consumo(N)*R/sum(Omega);
    Pag_Pacientes=Omega.*Aux;
end
    
Total=Pag_Impacientes+Pag_Pacientes;

% Autarquia
Impaciente_Autarquia=A*u(Dotacao_Economia/N);
Paciente_Autarquia=u(R*Dotacao_Economia/N);
BemEstarAutarquia=P*Impaciente_Autarquia+(1-P)*Paciente_Autarquia;

% Agora vamos calcular o bem estar no mecanismo ótimo calculado acima:

ConsumoOtimo=zeros(2^N,N); 
UtilidadeOtimo=zeros(2^N,1); 
OmegaOtimo=zeros(1,N);

BemEstarOtimo=0;

for e=1:2^N
    for j = 1:N
        OmegaOtimo(1,N-j+1) = floor(mod(e/(2^(j-1)),2));
    end
    % Aqui, pra cada omega calculado acima, colocamos o consumo ótimo e
    % assim por diante...
    ConsumoOtimo(e,:)=Consumo_Otimo(A, P, delta, R, OmegaOtimo, Dotacao_Economia, Lambda, F);
    UtilidadeOtimo(e,1)=(ones(1,N)-OmegaOtimo)*A*u(ConsumoOtimo(e,:))'+(OmegaOtimo)*u(ConsumoOtimo(e,:))';
    BemEstarOtimo=BemEstarOtimo+P^(N-sum(OmegaOtimo))*(1-P)^(sum(OmegaOtimo))*UtilidadeOtimo(e,1)/N;
end

BemEstarCorrida=UtilidadeOtimo(2^N,1)/N;

% Agora a mesma coisa para divesos estados:
clear all;
close all;
clc;

% Dados
A=1;
Ngrid=2:10;
R=1.05;
delta=2;
P=0.5;
a=1.11;
e=1;

% Utilidade
u=@(c)((c.^(1-delta))/(1-delta));

% Chute inicial para o multiplicador (seguimos o algoritmo proposto pelas
% notas):

for joao=1:length(Ngrid)

N=Ngrid(joao);    
Lambda=10;

dist=1;
tol=10^(-6);

itmax=10^(6);
it=0;

while it<itmax && dist>tol
    
      it=it+1;
      
      % Escrevendo gamma como no handout:
      
      gamma=abs(((1+Lambda/(1-P))/(A-Lambda/P))^(1/delta));
      
      M=N+1;
      F=zeros(M,M);
      
    % Obs:. Para não fazermos confusão com os índices, seguirei o padrão
    % estabelecido pelo handout. Desse modo utilizarei "i" para representar
    % o agente, e "j" para representar o anuncio.
    
    % Aqui vou calcular o F ótimo para o problema do agente após a
    % realização de um certo (w), conforme feito no handout.
    
    for j=1:M
        F(j,M)=gamma*(j-1)*R^(1/delta-1);
    end
    
    % Agora vou completar os temos da matriz F, a idéia e resolver o
    % problema retroativamente utilizando a informação calculada na parte
    % anterior.
      
      for k=1:M-1;
        j=M-k;
        for i=1:M-k
            F(i,j)=(P*((F(i,j+1)+1)^delta)+(1-P)*(F(i+1,j+1)^delta))^(1/delta);
        end    
      end
      
          
    % Agora preciso ver qual seria o ganho do agente de desviar para
    % podermos atualizar o multiplicador e iterar até chegar no ótimo:
      
      Ganho=zeros(M,M); 

    % Vamos fazer da mesma forma como fizemos para calcular F. Aqui o ganho
    % de desviar para o último indivíduo é equivalente ao que ele ganharia
    % em t=1 se desviasse do mecanismo calculado acima.
    
    for j=1:M 
        Ganho(j,M) = (P/(1-P))*(j-1)^delta*R^(1-delta);
    end
      
    % Completando a matriz:  
    
    % O calculo do ganho foi baseado novamente no handout.
      for k=1:M-1 
          j=M-k;
        for i=1:M-k
            Aux=Ganho(i,j+1)*(F(i,j+1)^(1-delta));
            if isnan(Aux)==1
                Aux=0;
            end    
            Ganho(i,j)=P*(Aux-1)*(1+F(i,j+1)^(delta-1))+(1-P)*(Ganho(i+1,j+1));
        end    
      end
      
    % Agora é só calcular o novo lambda, seguinda a metodologia da função
    % penalidade, verificar a convergência do memso e sair para o abraço
    
    lambda=exp(Ganho(1,1))*Lambda;
    dist=abs(Lambda-lambda);
    Lambda=lambda; 
end



% Autarquia
Dotacao_Economia=e*N;
Impaciente_Autarquia=A*u(Dotacao_Economia/N);
Paciente_Autarquia=u(R*Dotacao_Economia/N);
BemEstarAutarquia(joao)=P*Impaciente_Autarquia+(1-P)*Paciente_Autarquia;

% Agora vamos ver qual é o bem estar no mecanismo ótimo calculado acima:

ConsumoOtimo=zeros(2^N,N); 
UtilidadeOtimo=zeros(2^N,1); 
OmegaOtimo=zeros(1,N);

BemEstarOtimo(joao)=0;

for e=1:2^N
    for j = 1:N
        OmegaOtimo(1,N-j+1) = floor(mod(e/(2^(j-1)),2));
    end
    % Aqui pra cada omega calculado acima eu calculo o consumo ótimo e
    % assim por diante...
    ConsumoOtimo(e,:)=Consumo_Otimo(A, P, delta, R, OmegaOtimo, Dotacao_Economia, Lambda, F);
    UtilidadeOtimo(e,1)=(ones(1,N)-OmegaOtimo)*A*u(ConsumoOtimo(e,:))'+(OmegaOtimo)*u(ConsumoOtimo(e,:))';
    BemEstarOtimo(joao)=BemEstarOtimo(joao)+P^(N-sum(OmegaOtimo))*(1-P)^(sum(OmegaOtimo))*UtilidadeOtimo(e,1)/N;
end

BemEstarCorrida(joao)=UtilidadeOtimo(2^N,1)/N;
end

%% Plotando resultados
plot(Ngrid,BemEstarAutarquia,'red'); hold on;
plot(Ngrid,BemEstarOtimo, 'blue')
plot(Ngrid,BemEstarCorrida,'green')
hold off
legend('autarquia', 'mecanismo', 'bank-run')
ylabel('Bem-estar','FontSize',16)
xlabel('N','Fontsize',16)
title('Comparativo de utilidades', 'FontSize', 16)