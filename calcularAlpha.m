%% Cálculo do coeficiente de difusividade térmica (alpha)
% Alpha necessita ser recalculado em cada tempo j
% Input: Tj, temperatura média na seção, no tempo j
%        Te, temperatura de equilíbrio
%        T0, temperatura inicial
%        Xd, umidade média adimensional, na seção, no tempo j
%        j, identificador de tempo


function alpha = calcularAlpha(Tj, Te, T0, Xd)
%Temperatura média atual
% Td = (Tj - Te) / (T0 - Te); %Temperatura média adimensional, utilizada somente na eq. 21

% Difusividade Térmica (alpha) 
% Equação 20 e Tabela 2 Mariani et al. 2008.
A1 = 6.3197e-11 ;
A2 = 0.3027;
A3 = -0.3 ;
alpha = (A1 / (A2^Xd + A3));

% %Equação 21 e Tabela 3 Mariani et al. 2008.
% A1 = 5.4040e-11 ;
% A2 = 0.4758;
% A3 = -0.4741 ;
% A4 = 0.2465e10;
% alpha = (A1 / (A2^Xd + A3)) + A4/Td;
end