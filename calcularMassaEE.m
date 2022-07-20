%% Cálculo de massa
% Constroi a matriz de umidade por TDMA e o raio, considerando o
% encolhimento da banana durante aa secagem. 
% (OBS: Aplica refinamento na malha)
% Input: Xe,  umidade de equilíbrio
%        X0,  umidade incial
%        dt,  passo no tempo,
%        nt,  número de passos no tempo
%        R0,  raio inicial
%        nr,  número de nós na malha
%        Def, coeficiente difusivo
%        hm,  coeficiente convectivo de transferência de massa
%        f,   fator de refino da malha

% Output: X,  matriz de umidade
%         Xd, vetor de umidade média adimensional, na seção, no tempo j
%         R,  vetor de raios (totais), em cada do tempo j


function [X, Xd, R] = calcularMassaEE(Xe, X0, dt, nt, R0, nr, Def, hm, f)

X = zeros(f*nr,nt) ; % Inicialização da matriz
X(:,1) = X0; % Atribuição da umidade inicial no primeiro instante de tempo

%% Inicialização de variáveis

Xd = zeros(1,nt); %Registra o valor da umidade média adimensional da seção para cada tempo j
R = zeros(1, nt); %Registra tos os raios durante a secagem

%% Resolução
for j = 1:nt-1 % Para cada tempo j
    
    %Umidade média adimensional, na seção, no tempo j
    Xj = mean(X(:,j));
    Xd(j) = (Xj - Xe) / (X0 - Xe) ;
    
    %Cálculo do Raio
    R(j) = calcularRaio(R0, nr, Xj, Xe);
    r = linspace(0,R(j),f*nr);
    dr = R(j) / (f*nr);
    
    for i = f*nr:-1:1 % Para cada ponto i ao longo do raio
        
        % Atribuição de coeficientes
        if i == 1 %r = 0
            b = 2*Def*dt/(dr^2);
            X(i,j+1) = 2*b*X(i+1,j+1) + (1-2*b)*X(i,j);
            
        elseif i == f*nr %r = R
            b = Def/(hm*dr);
            X(i,j+1) = (Xe + b*X(i-1,j))/(1+b);
%             Xfin = X(i,j+1)
%             X2 = X(i, j)
            
            
        else %0 < r < R
            b = Def*dt/(r(i)*(dr^2));
            rp = (r(i+1) + r(i))/2;
            rm = (r(i-1) + r(i))/2;
            X(i,j+1) = b*rp*X(i+1,j+1) + (1-b*(rm+rp))*X(i,j) + b*rm*X(i-1,j);
        end
        
        %   plot((j-1)*dt/3600, X(1,j),'*')
        %   hold on
    end
end

%Umidade média adimensional, na seção, no tempo nt
Xj = mean(X(:,nt)) ;
Xd(nt) = (Xj - Xe) / (X0 - Xe) ;

%Cálculo do raio no instante final
R(nt) = (0.4721 + 0.1819*Xe + 0.1819*(Xj-Xe))*R0;
end
