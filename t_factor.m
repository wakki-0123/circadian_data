function t_factor(factor12000,t_kyoukai_plus,t_kyoukai_minus,t_atai)
figure;
semilogx(factor12000,t_kyoukai_plus,'-r');
xlim([0 13000]);
hold on;
semilogx(factor12000,t_kyoukai_minus,'-r');
xlim([0 13000]);
hold on;
semilogx(factor12000,t_atai,'-k');
xlim([0 13000]);

