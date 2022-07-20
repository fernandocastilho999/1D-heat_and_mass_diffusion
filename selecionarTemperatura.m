%% Selecionar Temperatura
% Constroi um vetor, a partir das temperaturas calculadas numéricamente,
% por meio de interpolação. Nesse vetor, cada temperatura refere-se a um 
% instante de tempo da medição da temperatura experimental, permitindo
% a comparação direta entre valores de temperatura numéricos e experiemntais.
%
% Input: T0,   temperatura incial
%        tnum, tempo de avaliação da temperatura numérica
%        texp, tempo de medição da temperatura experimental (tabelado)
%        nt,   número de passos no tempo
%        dt,   tamanho do passo no tempo
%
% Output: Tnum, vetor de temperaturas calculadas numéricamente, interpoladas
% para o tempo de medição experimental
%
% Método: Devido a malha empregada, o numero de instantes no tempo
% calculados numericamente (tnum) é maior do que o número de instantes no tempo de
% medição experimental (texp). Por esse motivo, essa função avalia qual tnum mais se
% aproxima de cada texp e, em seguida, faz uma interpolação linear entre a temperatura
% do tempo imediatamente anterior e a temperatura do tempo imediatamente posterior.

function Tnum = selecionarTemperatura(T0, tnum, texp, nt, dt)

% Inicialização do vetor
Tnum = zeros(1,135); % 135 texp, no total

cont = 1; % Variável de controle. Permite passar por cada um dos 135 texp. 

for j = 1:nt-1 % Avalia cada tnum(j), qual se aproxima mais do texp atual
    a = tnum(j)/3600; % tempo numérico atual, em h
    b = tnum(j+1)/3600; % próximo tempo numérico, em h
    c = T0(j); % temperatura no tempo atual
    d = T0(j+1); % temperatura no próximo tempo
    e = texp(cont); % tempo experimental avaliado, em h
    
    %interpoalação
    if (abs(e-a) <= dt/3600) % Se tnum(j) - texp é menor que o passo no tempo (h), interpola a temperatura correspondente
        Tnum(cont) = ((e-a)*(d-c)+c*(b-a))/(b-a);  % interpolação linear
        if cont < 135
            cont = cont + 1; % Passa a avliar o próximo texp da lista
        else
            break % Todos os 135 tempos experimentais já foram avaliados
        end
    end
end