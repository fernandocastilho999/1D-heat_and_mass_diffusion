%% Cálculo de massa
% Constroi a matriz de umidade por TDMA e o raio, considerando o
% encolhimento da banana durante aa secagem.
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


function [X, Xd, R] = calcularMassaDF(Xe, X0, dt, nt, R0, nr, Def, hm, f)

X = zeros(f*nr,nt) ; % Inicialização da matriz
X(:,1) = X0; % Atribuição da umidade inicial no primeiro instante de tempo

%% Inicialização de variáveis
% -a.X(i-1,j+1) + b.X(i,j+1) - c.X(i+1,j+1) = d
a = zeros(f*nr,1);
b = zeros(f*nr,1);
c = zeros(f*nr,1);
d = zeros(f*nr,1);
alfa = zeros(f*nr,1);
S = zeros(f*nr,1);

Xd = zeros(1,nt); %Registra o valor da umidade média adimensional da seção para cada tempo j
R = zeros(1, nt); %Registra tos os raios durante a secagem

%% Resolução

for j = 1:nt-1 % Para cada tempo j
    %Umidade média adimensional, na seção, no tempo j
    Xj = mean(X(:,j)) ;
    Xd(j) = (Xj - Xe) / (X0 - Xe) ;
    
    %Cálculo do raio para aquele instante de tempo devido ao encolhimento da banana
    R(j) = (0.4721 + 0.1819*Xe + 0.1819*(Xj-Xe))*R0;
    r = linspace(0,R(j),f*nr); % Vetor de raios igualmente espaçados, mantendo o número de nós constantes
    dr = R(j)/(f*nr); % Distância entre nós, na malha
    
    
    if j==1 % Calcula parqa j+1, ou seja, j=2, pelo método implícito, por TDMA
        for i = 1:f*nr % Para cada ponto i ao longo do raio
            
            % Atribuição de coeficientes
            if i == 1 %r = 0
                b(i) = 1 + 2*Def*dt/dr^2;
                c(i) = 2*Def*dt/dr^2;
                d(i) = (2*Def*dt/dr^2)*X(i+1,j) + (1-(2*Def*dt/dr^2))* X(i,j);
                alfa(i) = b(i);
                S(i) = d(i);
            elseif i == f*nr %r = R
                a(i) = Def/dr;
                b(i) = (Def/dr) + hm;
                d(i) = hm*Xe;
            else %0 < r < R
                bm = Def*dt/(2*r(i)*dr^2);
                rp = (r(i+1) + r(i))/2;
                rm = (r(i-1) + r(i))/2;
                a(i) = bm*rm;
                b(i) = (1 + bm*(rp+rm));
                c(i) = bm*rp;
                d(i) = bm*rp*X(i+1,j) + (1-bm*(rm+rp))*X(i,j) + bm*rm*X(i-1,j);
            end
            
            %Solução da matriz
            if i~=1
                alfa(i)=b(i)-a(i)*c(i-1)/alfa(i-1);
                S(i)=d(i)+a(i)*S(i-1)/alfa(i-1);
            end
        end
        
        %Solução do último nó
        X(f*nr,j+1)=S(f*nr)/alfa(f*nr);
        
        %Substituição regressiva
        for i=f*nr-1:-1:1
            X(i,j+1)=(S(i)+c(i)*X(i+1,j+1))/alfa(i);
        end
        
    else % Para todos os j restantes, ou seja, j >= 3, calcula por Dufort-Frankel
        
        for i = 1:f*nr % Para cada ponto da malha
            
            if i ==1 % r=0
                b= 2*Def*dt/dr^2;
                X(i,j+1) = ((1-b)*X(i,j-1) + 2*b*X(i+1,j))/(1+b);
                
            elseif i == f*nr % r=R
                b =Def/dr;
                X(i,j+1) = (b*X(i-1,j+1) +hm*Xe)/(b+hm);
                
            else % 0 < r < R
                b = 2*Def*dt/ (r(i) * dr^2)  ; % constante para facilitar os cálculos
                rm = (r(i-1) + r(i))/2 ; % rm = r[-1/2]
                rp = (r(i) + r(i+1))/2 ; % rp = r[+1/2]
                X(i,j+1) = (b*rm*X(i-1,j)+(1-b*r(i))*X(i,j-1)+b*rp*X(i+1,j))/(1+b*r(i));
                
            end
        end
%         plot((j-1)*dt, X(1,j),'*')
%         hold on
    end
end
end
