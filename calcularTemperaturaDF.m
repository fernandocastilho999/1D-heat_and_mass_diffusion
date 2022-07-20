%% Temperatura, coeficientes e propriedades que variam com tempo e raio
% Constroi a matriz de temperatura por TDMA. Chama as funções de cálculo de propriedades, que dependendem do
% tempo e do raio.
% Input: Te,  temperatura de equilíbrio
%        T0,  temperatura incial
%        X,   matriz de umidade
%        Xd,  vetor de umidade média adimensional, na seção, em cada tempo j
%        dt,  passo no tempo,
%        nt,  número de passos no tempo
%        R,   vetor de raios (totais), em cada do tempo j
%        nr,  número de nós na malha
%        v,   velocidade do vento
%        hfg, calor latente de vaporização da água
%        cv,  calor específico do vapor de água
%        rhos, densidade do sólido


function [T, alpha] = calcularTemperaturaDF(Te, T0, X, Xd, dt, nt, R, nr, v, hfg, cv, rhos)

T = zeros(nr,nt); % Inicialização da matriz de temperatura
T(:,1) = T0; % Atribuição da temperatura inicial no primeiro instante de tempo
alpha = zeros(1, nt); % Inicialização do vetor de coeficiente de difusividade térmica

%% Inicialização de variáveis
% -a.T(i-1,j+1) + b.T(i,j+1) - c.T(i+1,j+1) = d
a = zeros(nr,1);
b = zeros(nr,1);
c = zeros(nr,1);
d = zeros(nr,1);
beta = zeros(nr,1);
S = zeros(nr,1);
% dxdt = zeros(1,nt);
% h = zeros(1,nt);

%% Resolução
for j = 1:nt-1
    
    %Cálculo do raio para aquele instante de tempo devido ao encolhimento da banana
    r = linspace(0,R(j),nr); % Vetor de raios igualmente espaçados, mantendo o número de nós constantes
    dr = R(j)/nr; % Distância entre nós, na malha
    
    %Cálculo do coeficiente de difusividade térmica (alpha) em cada instante j
    alpha(j) = calcularAlpha(mean(T(:,j)), Te, T0, Xd(j));
    
    Tj = 0;
    for i = 1:nr
        Tj = Tj + T(i,j)*r(i)/R(j);
    end
    
    if j == 1 % Calcula j+1, ou seja, j = 2 com o método implícito, por TDMA
        for i = 1:nr %Para cada nó i ao longo do raio
            
            % Atribuição de coeficientes
            if i == 1 %r = 0
                
                b(i) = 1 + (2*dt*alpha(j)/dr^2);
                c(i) = 2*dt*alpha(j)/dr^2;
                d(i) = (2*alpha(j)*dt/dr^2)*T(i+1,j) + (1-(2*alpha(j)*dt/dr^2))* T(i,j);
                beta(i) = b(i);
                S(i) = d(i);
                
            elseif i == nr %r = R
                
                %Propriedades
                [h, k] = calcularPropriedades(T(i,j), mean(X(:,j)), v, alpha(j), R(j));
                dxdt = (mean(X(:,j+1)) - mean(X(:,j)))/dt;
                
                a(i) = 1;
                b(i) = 1 + h*dr/k + (rhos*(dr^2)/k)*dxdt*cv;
                d(i) = h*dr*Te/k - (rhos*(dr^2)/k)*dxdt*(hfg-cv*Te);
                
            else % 0 < r < R
                
                bm = alpha(j)*dt/(2*r(i)*dr^2);
                rp = (r(i+1) + r(i))/2;
                rm = (r(i-1) + r(i))/2;
                
                a(i) = bm*rm;
                b(i) = 1 + bm*(rp+rm);
                c(i) = bm*rp;
                d(i) = bm*rp*T(i+1,j) + (1-bm*(rm+rp))*T(i,j) + bm*rm*T(i-1,j);
                
            end
            
            if i~=1
                beta(i) = b(i)-a(i)*c(i-1)/beta(i-1);
                S(i)= d(i)+a(i)*S(i-1)/beta(i-1);
            end
        end
        
        T(nr,j+1)=S(nr)/beta(nr); %Solução do último nó
        
        % Substituição retroativa
        for i=nr-1:-1:1
            T(i,j+1)=(S(i)+c(i)*T(i+1,j+1))/beta(i);
        end
        
    else  % Para todos os j restantes, ou seja, j >= 3, calcula por Dufort-Frankel
        
        for i = 1:nr % Para cada ponto da malha
            
            if i == 1 % r=0
                
                b = 2*alpha(j)*dt/dr^2;
                T(i,j+1) = ((1-b)*T(i,j-1)+2*b*T(i+1,j))/(1+b);
                
            elseif i == nr % r=R
                
                %Propriedades
                [h, k] = calcularPropriedades(T(i,j), Tj, v, alpha(j), R(j));
                dxdt = (mean(X(:,j+1)) - mean(X(:,j)))/dt;
                
                b = 1 + dr*h/k + (rhos*(dr^2)/k)*dxdt*cv;
                c = dr*h*Te/k - (rhos*(dr^2)/k)*dxdt*(hfg-cv*Te);
                T(i,j+1) = (T(i-1,j+1) + c)/b;
                
            else % 0 < r < R
                
                b = 2*alpha(j)*dt/ (r(i) * dr^2 );
                rm = (r(i-1) + r(i))/2;
                rp = (r(i) + r(i+1))/2;
                T(i,j+1) = (b*rp*T(i+1,j) + (1-b*r(i))*T(i,j-1) +b*rm*T(i-1,j))/(1+b*r(i));
                
            end
        end
        %     plot((j-1)*dt, T(1,j),'*')
        %     hold on
    end
    
end
end

