function t_factor(t_atai)
    factor = 9993;
    num = 7;

   % t_kyoukai_plus = 4.057735 * ones(1, 10000);  
   % t_kyoukai_minus = -4.057735 * ones(1, 10000); 
% Fuzzy
    t_kyoukai_plus = 2.433081 * ones(1, factor);  
    t_kyoukai_minus = -2.433081 * ones(1, factor); 

   % MSE
    % t_kyoukai_plus = 3.195322 * ones(1, 1000);  
    % t_kyoukai_minus = -3.195322 * ones(1, 1000); 
    
time_length = factor * 5; % 全部の区間の秒数

% 時間スケールを計算
time_s = zeros(1, factor);
time = zeros(1, factor);
for i = 1:factor
    time_s(i) = factor / i; % 合計サンプルの個数
    time(i) = (time_length / time_s(i)) + (5 * num); % タイムスケール
end


    figure;
    plot(time, t_atai,'Color',[0 0 0],'LineWidth',4);      % t_atai をプロット
    hold on;
    plot(time, t_kyoukai_plus,'--r','LineWidth',3);   % t_kyoukai_plus をプロット
    plot(time, t_kyoukai_minus,'--b','LineWidth',3);  % t_kyoukai_minus をプロット
    hold off;
    set(gca, 'XScale', 'log');   % X軸を対数スケールに設定
    grid on
    xlim([time(1) time(factor)])
    xlabel('Time Scale [sec]');             % X軸ラベルを設定
    ylabel('t-value');             % Y軸ラベルを設定
    lgd = legend('t-value', 'Upper threshold  t-value corresponding to q<0.05', 'lower threshold  t-value corresponding to q<0.05','Location', 'southeast');  % 凡例を追加
    lgd.FontSize = 18;
    ax = gca;
    ax.FontSize = 45;
    %title('Plot of t\_atai, t\_kyoukai\_plus, and t\_kyoukai\_minus');  % タイトルを追加
end
