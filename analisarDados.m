%% Análise de dados
% Calcula o MAE, o MAPE, o RMSE e o R² para a curva de temperatura.
% Input: Tr0,  vetor de temperaturas no centro da banana (r=0 ou i=1)
%        Texp, vetor de temperaturas experimentais (Pérez, 1998)
%
% Output: MAE, erro absoluto médio
%         MAPE, erro percentual absoluto médio
%         RMSE, raiz quadrada do erro médio
%         R2, coeficiente de Pearson

function [MAE, MAPE, RMSE, R2] = analisarDados(Tr0, Texp)

% Inicialização de variáveis
MAE = 0; % MAE = SOMA[j=1:nt](abs(Texp(j) - Tr0(j)))/nt
MAPE = 0; % MAPE = SOMA[j=1:nt](abs((Texp(j) - Tr0(j))/Texp(j)))*100/nt
RMSE = 0; % RMSE = SQRT(SOMA[j=1:nt]((Texp(j) - Tr0(j))^2)/nt)
%R2 = 0; % R2 = 1 - SOMA[j=1:nt](Texp(j) - Tr0(j))^2 / SOMA[j=1:nt](Texp(j) - mean(Texp(j)))^2

% Avalia cada um dos 135 valores experimentais
for j = 1:135 
    MAE = MAE + abs(Texp(j) - Tr0(j));
    MAPE = MAPE + abs((Texp(j) - Tr0(j))/Texp(j));
    RMSE = RMSE + (Texp(j) - Tr0(j))^2;
end

MAE = MAE / 135;
MAPE = MAPE * 100/135;
RMSE = sqrt(RMSE / 135);
R2 = corrcoef(Tr0, Texp);

end