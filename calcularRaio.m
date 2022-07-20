%% Cálculo do Raio com encolhimento
% 


function [R, r, dr] = calcularRaio(R0, nr, Xj, Xe)  
    %Cálculo do raio para naquele instante de tempo devido ao encolhimento da banana
    R = (0.4721 + 0.1819*Xe + 0.1819*(Xj-Xe))*R0;
    r = linspace(0,R,nr); % Vetor de raios igualmente espaçados, mantendo o número de nós constantes
    dr = R/nr; % Distância entre nós, na malha
end