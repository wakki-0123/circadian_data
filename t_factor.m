function t_factor(t_atai)
   % t_kyoukai_plus = 4.057735 * ones(1, 10000);  
   % t_kyoukai_minus = -4.057735 * ones(1, 10000); 
% Fuzzy
    t_kyoukai_plus = 4.614038 * ones(1, 10000);  
    t_kyoukai_minus = -4.614038 * ones(1, 10000); 

   % MSE
    % t_kyoukai_plus = 3.195322 * ones(1, 1000);  
    % t_kyoukai_minus = -3.195322 * ones(1, 1000); 
    
time_length = 10000 * 5; % 全部の区間の秒数
factor = 10000;
% 時間スケールを計算
time_s = zeros(1, factor);
time = zeros(1, factor);
for i = 1:factor
    time_s(i) = 10000 / i; % 合計サンプルの個数
    time(i) = time_length / time_s(i); % タイムスケール
end


    figure;
    plot(time, t_atai,'Color',[0 0 0]);      % t_atai をプロット
    hold on;
    plot(time, t_kyoukai_plus,'--r');   % t_kyoukai_plus をプロット
    plot(time, t_kyoukai_minus,'--b');  % t_kyoukai_minus をプロット
    hold off;
    set(gca, 'XScale', 'log');   % X軸を対数スケールに設定
    grid on
    xlim([time(1) time(10000)])
    xlabel('Time Scale [sec]');             % X軸ラベルを設定
    ylabel('t-Value');             % Y軸ラベルを設定
    lgd = legend('t-value', 'Upper threshold  t-Value corresponding to q<0.05', 'lower threshold  t-Value corresponding to q<0.05','Location', 'southeast');  % 凡例を追加
    lgd.FontSize = 27;
    ax = gca;
    ax.FontSize = 45;
    %title('Plot of t\_atai, t\_kyoukai\_plus, and t\_kyoukai\_minus');  % タイトルを追加
end
