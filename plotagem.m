function plotagem(T, Xd, alpha, temp, nt, experimental)

% Gráfico comparativo entre dados calculados numericamente e dados experimentais tabelelados por Pérez, 1998
subplot(1,3,1)
plot( temp(:)/3600, T(1,:), 'k-', experimental(:,1), experimental(:,2), 'rx')   
xlim([0 140]) ; ylim([15 40]) ; 
grid on
title (["Temperatura Calculada x Temperatura Experimental"]);
xlabel ("Tempo (h)");
ylabel ("Temperatura no centro R0 - T");
hold on

% Curva de secagem, umidade adimensional, média na seção x tempo de secagem
subplot(1,3,2)
plot( temp(1:nt-1)/3600, Xd(1:nt-1), 'k-')
xlim([0 140]) ; ylim([0 1]) ; 
grid on
title (["Curva de secagem"]);
xlabel ("Tempo (h)");
ylabel ("Umidade média adimensional na seção");
hold on

% Curva de secagem, umidade adimensional, média na seção x ln(alpha)
subplot(1,3,3)
plot( Xd(1:nt-1), log(alpha(1:nt-1)), 'k-')
xlim([0.82 1]) ; ylim([-20.5 -15.5]) ;
grid on
title (["Curva de secagem  por ln(α)"]);
xlabel ("Umidade média adimensional na seção");
ylabel ("ln(α)");
hold on


end