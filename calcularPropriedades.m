%% Propriedades que variam com tempo e raio
% Essas prorpiedades necessitam ser recalculadas em cada tempo j e nó i
% Input:   Tij, temperatura no nó i, no tempo j
%           Xj, umidade média na seção, no tempo j
%            v, velocidade do ar
%        alpha, coefciente de difusividade no tempo j
%            r, é o raio da banana no tempo j
%
% Output:    h, coeficiente convectivo
%            k, condutividade térmica do ar


function [h, k] = calcularPropriedades(Tij, Xj, v, alpha, rj)

%Temperatura em Kelvin
Tk = Tij+ 273.15 ; 

%Densidade do ar
P = 101325 ;
R = 287.058 ;
rho = P / R / Tk ; %densidade do ar seco

UA = 0.357 * 30.4e-3; %Umidade absoluta kg/m³
% UR = 35,7%, de acordo com a dissertação do Perez

rhou = rho + UA; %Densidade do ar úmido 

%Viscosidade do ar
muo = 1.7894e-5 ;
To = 273.11 ;
S = 110.15 ;
mu = muo * (Tk/To)^(3/2) * (To+S)/(Tk+S) ;

%Umidade em base úmida (considerando o peso específico da água = 1000 kg/m³)
Xjw = Xj / (1 + Xj);

%Condutividade térmica da banana (?)
k = 0.148 + 0.493*Xjw;

% %Condutividade térmica do ar
ka = 0.02588;

%Número de Prandtl
% Pr = mu/(rhou*alpha);
Pr = mu/(rhou*2.208e-5);

%Número de Reynolds
Re = rhou*v*2*rj/mu ;

%Número de Nusselt
Nu = 0.97 + 0.68* Re^0.52 * Pr^(1/3) ;

%Coeficiente convectivo
h = ka*Nu/(2*rj) ;

end