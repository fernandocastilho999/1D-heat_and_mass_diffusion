clear all
close all
clc

experimental = load("Temp_EXP.mat").varitempexp; %Importa os dados experimentais de Pérez, 1998

%Parâmetros - Caso 1
Te = 29.9; % Temperatura de equilíbrio
hm = 3.19e-7; % Coeficiente convectivo de transferência de massa
Def = 1.48e-10; % Coeficiente difusivo
Xe = 0.1428; % Umidade de equilíbrio
v = 0.38;  % Velocidade do ar
hfg = 2430750; % Calor latente de vaporização da água
cv = 1902.3; % Calor específico do vapor de água
t = 121.9*3600; % Tempo total
rhos = 1970; % Densidade do sólido

% Condições Iniciais
T0 = 19.1;
X0 = 3.43;
R0 = 0.015;
tf = t;

%Definição da malha
% nr = 50;
% nt = 200*3600;
nr = 8;
nt = 500000;
dt = tf/nt;

%Inicialização das variáveis
temp = linspace(0,tf,nt); % Malha tempo

%% Cálculo de massa 
%Umidade X, Umidade média adimensional Xd e Raio R em cada tempo j

f = 50; %Fator de refino da malha espacial

[X_CN, Xd_CN, R_CN] = calcularMassaCN(Xe, X0, dt, nt, R0, nr, Def, hm, f); %Crank-Nicolson
[X_EE, Xd_EE, R_EE] = calcularMassaEE(Xe, X0, dt, nt, R0, nr, Def, hm, f); % Euler Explícito
[X_DF, Xd_DF, R_DF] = calcularMassaDF(Xe, X0, dt, nt, R0, nr, Def, hm, f); % Dufort-Frankel

%% Cálculo de temperatura e difusividade térmica
%Temperatura T e coeficiente de difusividade térmica alpha

[T_CN, alpha_CN] = calcularTemperaturaCN(Te, T0, X_CN, Xd_CN, dt, nt, R_CN, nr, v, hfg, cv, rhos); % Crnak-Nicolson 
[T_EE, alpha_EE] = calcularTemperaturaEE(Te, T0, X_EE, Xd_EE, dt, nt, R_EE, nr, v, hfg, cv, rhos); % Euler Explícito 
[T_DF, alpha_DF] = calcularTemperaturaDF(Te, T0, X_DF, Xd_DF, dt, nt, R_DF, nr, v, hfg, cv, rhos); % Dufort-Frankel 

%% Analise estatística dos dados obtidos

% Seleciona as temperaturas correspondentes aos dados experimentais
Tnum_CN = selecionarTemperatura(T_CN(1,:), temp, experimental(:,1), nt, dt); 
Tnum_EE = selecionarTemperatura(T_EE(1,:), temp, experimental(:,1), nt, dt); 
Tnum_DF = selecionarTemperatura(T_DF(1,:), temp, experimental(:,1), nt, dt); 

[MAE_CN, MAPE_CN, RMSE_CN, R2_CN] = analisarDados(Tnum_CN, experimental(:,2));
[MAE_EE, MAPE_EE, RMSE_EE, R2_EE] = analisarDados(Tnum_EE, experimental(:,2));
[MAE_DF, MAPE_DF, RMSE_DF, R2_DF] = analisarDados(Tnum_DF, experimental(:,2));

%% Plotagem dos gráficos

plotagem3(T_CN, Xd_CN, alpha_CN, T_EE, Xd_EE, alpha_EE, T_DF, Xd_DF, alpha_DF, temp, nt, experimental);

% plotagem (T_CN, Xd_CN, alpha_CN, temp, nt, experimental);
% plotagem (T_EE, Xd_EE, alpha_EE, temp, nt, experimental);
% plotagem (T_DF, Xd_DF, alpha_DF, temp, nt, experimental);
