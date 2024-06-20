function sleep_t_factor(sleep_correlation_coefficients,scale)
   % t_kyoukai_plus = 4.057735 * ones(1, 10000);  
   % t_kyoukai_minus = -4.057735 * ones(1, 10000); 
% Fuzzy
    t_kyoukai_plus = 0.818182 * ones(1, scale);  
    t_kyoukai_minus = -0.818182 * ones(1, scale); 

   % MSE
    % t_kyoukai_plus = 3.195322 * ones(1, 1000);  
    % t_kyoukai_minus = -3.195322 * ones(1, 1000); 
    


time_length = scale * 5; % 全部の区間の秒数
%factor = 10000;
% 時間スケールを計算
time_s = zeros(1, scale);
time = zeros(1, scale);
for i = 1:scale
    time_s(i) = scale / i; % 合計サンプルの個数
    time(i) = time_length / time_s(i); % タイムスケール
end

    figure;
    plot(time,sleep_correlation_coefficients,'Color',[0 0 0],'LineWidth',5);      % t_atai をプロット
    hold on;
    plot(time,t_kyoukai_plus,'--r','LineWidth',3);   % t_kyoukai_plus をプロット
    plot(time,t_kyoukai_minus,'--b','LineWidth',3);  % t_kyoukai_minus をプロット
    set(gca, 'XScale', 'log');
    xlim([time(1) time(scale)])
    hold off;
   
    grid on
    
    xlabel('Time Scale');             % X軸ラベルを設定
    ylabel('correlation coefficients');             % Y軸ラベルを設定
    lgd = legend('t-value', 'Upper threshold  t-value corresponding to q<0.05', 'lower threshold  t-value corresponding to q<0.05','Location', 'southeast');  % 凡例を追加
    lgd.FontSize = 20;
    ax = gca;
    ax.FontSize = 40;
    %title('Plot of t\_atai, t\_kyoukai\_plus, and t\_kyoukai\_minus');  % タイトルを追加
end
