%% Temperatura, coeficientes e propriedades que variam com tempo e raio
% Constroi a matriz de temperatura por TDMA. Chama as funções de cálculo de propriedades, que dependendem do
% tempo e do raio. GAUSS-SEIDEL
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


function [T, alpha] = calcularTemperaturaEE(Te, T0, X, Xd, dt, nt, R, nr, v, hfg, cv, rhos)

%% Inicialização de variáveis
T = zeros(nr,nt); % Inicialização da matriz de temperatura
T(:,1) = T0; % Atribuição da temperatura inicial no primeiro instante de tempo
alpha = zeros(1, nt); % Inicialização do vetor de coeficiente de difusividade térmica

%% Resolução
for j = 1:nt-1 %Para cada tempo j
    
    %Cálculo do raio para naquele instante de tempo devido ao encolhimento da banana
    r = linspace(0,R(j),nr); % Vetor de raios igualmente espaçados, mantendo o número de nós constantes
    dr = R(j)/nr; % Distância entre nós, na malha
    
    %Cálculo do coeficiente de difusividade térmica (alpha) em cada instante j
    alpha(j) = calcularAlpha(mean(T(:,j)), Te, T0, Xd(j));
    
    for i = nr:-1:1 %Para cada nó i ao longo do raio, começando de r = R e seguindo até r = 0
        
        if i == 1 %r = 0
            b = alpha(j)*dt/dr^2;
            T(i,j+1) = 2*b*T(i+1,j+1) + (1-2*b)*T(i,j);
            
        elseif i == nr %r = R
            % Propriedades necessárias: Difusividade términa (alpha), coeficiente convectivo (h) e conductividade térmica do ar(k)
            [h, k] = calcularPropriedades(T(i,j), mean(X(:,j)), v, alpha(j), R(j));
            dxdt = (mean(X(:,j+1)) - mean(X(:,j)))/dt;
            
            b = dr*h*Te/k - (rhos*(dr^2)/k)*dxdt*(hfg-cv*Te);
            cj = 1 + dr*h/k + (rhos*(dr^2)/k)*dxdt*cv;
            T(i,j+1) =  (T(i-1,j)+ b)/cj;
            
        else % 0 < r < R
            b = alpha(j)*dt/(r(i)*dr^2);
            rp = (r(i+1) + r(i))/2;
            rm = (r(i-1) + r(i))/2;
            T(i,j+1) = b*rp*T(i+1,j+1) + (1-b*(rp+rm))*T(i,j) + b*rm*T(i-1,j);
        end
        
    end
    %         plot((j-1)*dt/3600, T(1,j),'*')
    %         hold on
end
end