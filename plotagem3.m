function plotagem3(T_CN, Xd_CN, alpha_CN, T_EE, Xd_EE, alpha_EE, T_DF, Xd_DF, alpha_DF, temp, nt, experimental)

% Gráfico comparativo entre dados calculados numericamente e dados experimentais tabelelados por Pérez, 1998
subplot(1,3,1)
plot(experimental(:,1), experimental(:,2), 'kx', temp(:)/3600, T_CN(1,:), 'r-', temp(:)/3600, T_EE(1,:), 'b-', temp(:)/3600, T_DF(1,:), 'g-')   
xlim([0 140]) ; ylim([15 40]) ; 
grid on
title (["Temperatura Calculada x Temperatura Experimental"]);
xlabel ("Tempo (h)");
ylabel ("Temperatura no centro R0 - T");
legend('Experimental','CN', 'EE', 'DF');
hold on

% Curva de secagem, umidade adimensional, média na seção x tempo de secagem
subplot(1,3,2)
plot( temp(1:nt-1)/3600, Xd_CN(1:nt-1), 'r-', temp(1:nt-1)/3600, Xd_EE(1:nt-1), 'b-', temp(1:nt-1)/3600, Xd_DF(1:nt-1), 'g-')
xlim([0 140]) ; ylim([0 1]) ; 
grid on
title (["Curva de secagem"]);
xlabel ("Tempo (h)");
ylabel ("Umidade média adimensional na seção");
legend('CN', 'EE', 'DF');
hold on

% Curva de secagem, umidade adimensional, média na seção x ln(alpha)
subplot(1,3,3)
plot( Xd_CN(1:nt-1), log(alpha_CN(1:nt-1)), 'r-', Xd_EE(1:nt-1), log(alpha_EE(1:nt-1)), 'b-', Xd_DF(1:nt-1), log(alpha_DF(1:nt-1)), 'g-')
xlim([0.82 1]) ; ylim([-20.5 -15.5]) ;
grid on
title (["Curva de secagem  por ln(α)"]);
xlabel ("Umidade média adimensional na seção");
ylabel ("ln(α)");
legend('CN', 'EE', 'DF');
hold on


end