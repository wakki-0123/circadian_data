function t_factor(t_atai)
    t_kyoukai_plus = 4.057735 * ones(1, 10000);  
    t_kyoukai_minus = -4.057735 * ones(1, 10000); 
    figure;
    plot(1:10000, t_atai);      % t_atai をプロット
    hold on;
    plot(1:10000, t_kyoukai_plus,'-r');   % t_kyoukai_plus をプロット
    plot(1:10000, t_kyoukai_minus,'-b');  % t_kyoukai_minus をプロット
    hold off;
    set(gca, 'XScale', 'log');   % X軸を対数スケールに設定
    xlim([1 1000])
    xlabel('Index');             % X軸ラベルを設定
    ylabel('Value');             % Y軸ラベルを設定
    legend('t\_atai', 't\_kyoukai\_plus', 't\_kyoukai\_minus');  % 凡例を追加
    title('Plot of t\_atai, t\_kyoukai\_plus, and t\_kyoukai\_minus');  % タイトルを追加
end
