% function [p_values, correlation_coefficients, t_values] = sleep_scatter_0823(e)
%     % 被験者のテストスコア
%    %test_score = [0.886324022, 0.866727644, 0.911264649, 0.877716281, 0.885269186, 0.855124121, 0.871734582, 0.859731051, 0.900972038, 0.854969546];
%     test_score = [0.886324022, 0.866727644, 0.911264649, 0.877716281, 0.885269186, 0.855124121, 0.871734582, 0.859731051, 0.900972038, 0.854969546,0.8519555, 0.894203152];
%       % , 0.8519555, 0.894203152
%     entropy = e; % 12x10000 行列
%     num = 12; % 被験者の数
%     scale =1000;
% 
%     % 相関係数、p値、t値の初期化
%     correlation_coefficients = zeros(1, scale);
%     p_values = zeros(1, scale);
%     t_values = zeros(1, scale);
% 
%     % 各時間スケールのデータを処理
%     for i = 1:scale
%         % 各時間スケールごとにエントロピー値を取り出す
%         E_person = entropy(:, i); % 12x1 ベクトル
% 
%         % Spearmanの順位相関係数を計算
%         [correlation_coefficients(i), p_values(i)] = corr(E_person, test_score', 'Type', 'Spearman');
% 
%         % t値の計算
%         n = num; % サンプル数
%         t_values(i) = correlation_coefficients(i) * sqrt((n - 2) / (1 - correlation_coefficients(i)^2));
%     end
% 
%     % 相関係数の図の作成
%     figure;
%     scatter(1:scale, correlation_coefficients, 'filled');
%     xlabel('Time Scale');
%     ylabel('Spearman Correlation Coefficient');
%     title('Correlation Coefficients over Different Time Scales');
%     grid on;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %1121使用
% function [p_values, correlation_coefficients, t_values] = sleep_scatter_0823(e)
%     % 被験者のテストスコア
%     % test_score = [0.886324022, 0.866727644, 0.911264649, 0.877716281, 0.885269186, ...
%     %               0.855124121, 0.871734582, 0.859731051, 0.900972038, 0.854969546, ...
%     %               0.8519555, 0.894203152,0.905361808]; 
% % 15人
%     %test_score = [0.876847548,0.867301668,0.888202161,0.87524095,0.85704005,0.911264649,0.855918089,0.901118598,0.859571395,0.886905742,0.891167085,0.8519555,0.890027395,0.890409196,0.916248705,0.872267,0.870477,0.884756];
% 
%     test_score = [0.876847548,0.867301668,0.888202161,0.87524095,0.85704005,0.911264649,0.855918089,0.901118598,0.859571395,0.886905742,0.891167085,0.8519555,0.890027395,0.916248705,0.872267,0.884756];
%     % 14行目と17行目を除くための行インデックスを作成します
%     rows_to_exclude = [14, 17];
% 
% % 行インデックスを利用して、新しい配列を作成します
%     e1 = e(setdiff(1:size(e, 1), rows_to_exclude), :);  
% 
% 
%     % 1120現在
%     entropy = e1; % 12x10000 行列
%     num = 16; % 被験者の数
%     %scale = 9993; % スケール数を変更
%     scale = 600;
%     scale_bottom = 1 ;% スケールの底
% 
%     % 相関係数、p値、t値の初期化
%     correlation_coefficients = zeros(1, scale-scale_bottom+1);
%     p_values = zeros(1, scale-scale_bottom+1);
%     t_values = zeros(1, scale-scale_bottom+1);
%     num_2 = 7; % 省略したfactorの数(0926 現在は7)
% 
%     % 各時間スケールのデータを処理
%     for i = scale_bottom:1:scale
%         % 各時間スケールごとにエントロピー値を取り出す
%         E_person = entropy(1:num, i); % numx1 ベクトル
% 
%         % Spearmanの順位相関係数を計算
%         [correlation_coefficients(i - scale_bottom + 1), p_values(i - scale_bottom + 1)] = corr(E_person, test_score(1:num)', 'Type','Spearman');
% 
%         % t値の計算
%         n = num; % サンプル数
%         t_values(i - scale_bottom + 1) = correlation_coefficients(i - scale_bottom + 1) * sqrt((n - 2) / (1 - correlation_coefficients(i - scale_bottom + 1)^2));
%     end
% 
%     disp(size(correlation_coefficients))
% 
%     % 相関係数の図の作成
%     figure;
%     scatter((scale_bottom+num_2) * 5:5: (scale+num_2) * 5, correlation_coefficients, 'filled');
%     xlabel('Time Scale [s]');
%     ylabel('Spearman Correlation Coefficient');
%     title('Correlation Coefficients over Different Time Scales');
%     grid on;
% end

%%%%%%%%%%%%%%%%%%%

% 1210時点最新
function [p_values, correlation_coefficients, t_values, l, entropy] = sleep_scatter_0823(e)
    % 被験者のテストスコア 16人分,14,17の人を除いた
    test_score = [0.876847548,0.867301668,0.888202161,0.87524095,0.85704005,0.911264649,0.855918089,0.901118598,0.859571395,0.886905742,0.891167085,0.8519555,0.890027395,0.916248705,0.872267,0.884756];
    % 14行目と17行目を除くための行インデックスを作成します
    rows_to_exclude = [14, 17];

    % 行インデックスを利用して、新しい配列を作成します
    e1 = e(setdiff(1:size(e, 1), rows_to_exclude), :);

    entropy = e1; % 16x10000 行列
    %entropy = e;
    num = 16; % 被験者の数 元は18
    scale = 800; % スケール数
    scale_bottom = 1; % スケールの底

    % 相関係数、p値、t値の初期化
    correlation_coefficients = zeros(1, scale-scale_bottom+1);
    p_values = zeros(1, scale-scale_bottom+1);
    t_values = zeros(1, scale-scale_bottom+1);
    num_2 = 7; % 省略したfactorの数(0926 現在は7)

    % 各時間スケールのデータを処理
    for i = scale_bottom:1:scale
        % 各時間スケールごとにエントロピー値を取り出す
        E_person = entropy(1:num, i); % numx1 ベクトル

        % Spearmanの順位相関係数を計算
        [correlation_coefficients(i - scale_bottom + 1), p_values(i - scale_bottom + 1)] = corr(E_person, test_score(1:num)', 'Type','Spearman'); %Spearman

        % t値の計算
        n = num; % サンプル数
        t_values(i - scale_bottom + 1) = correlation_coefficients(i - scale_bottom + 1) * sqrt((n - 2) / (1 - correlation_coefficients(i - scale_bottom + 1)^2));
    end

    disp('point number')
    disp(size(correlation_coefficients));

    % p値が0.05未満のスケールを収集
    l = [];
    s = 0;
    for i = 1:length(p_values)
        if p_values(i) < 0.05
            s = s + 1;
            l(s) = i;
        end
    end

    % 相関係数の図の作成
    %%%%%%%%%%%%%
    % % scale time
    figure;
    hold on;
    time_scales = (scale_bottom+num_2) * 5:5: (scale+num_2) * 5;
    significant = p_values < 0.05;
    non_significant = ~significant;

    scatter(time_scales(significant), correlation_coefficients(significant), 'filled', 'b'); %r
    scatter(time_scales(non_significant), correlation_coefficients(non_significant), 'filled', 'b');

    % 境界線の追加
    yline(min(correlation_coefficients(significant)), '--r', 'p = 0.05', 'LabelHorizontalAlignment', 'left','LineWidth',5, 'FontSize', 25);

    xlabel('Time Scale [s]', 'FontSize', 25);
    ylabel('Spearman Correlation Coefficient', 'FontSize', 25);

    %title('Correlation Coefficients over Different Time Scales', 'FontSize', 18);
    %legend('p < 0.05', 'p \geq 0.05', 'FontSize', 14);

    ax = gca;
    ax.FontSize = 20;

    grid on;
    hold off;

    % 相関係数が最大および最小のタイムスケールを求める
    [max_corr, max_idx] = max(correlation_coefficients);
    [min_corr, min_idx] = min(correlation_coefficients);
    
    % 対応するエントロピー値を取得
    E_max = entropy(:, max_idx); % 16x1 ベクトル
    E_min = entropy(:, min_idx); % 16x1 ベクトル

    % 散布図を作成
    figure;
    scatter(E_max, test_score, 'filled', 'b');
    xlabel('Fuzzy Entropy', 'FontSize', 25);
    ylabel('Sleep Efficiency', 'FontSize', 25);
    title(['Scatter Plot at Max Correlation Time Scale ', num2str((max_idx + num_2) * 5)], 'FontSize', 25);
    ax = gca;
    ax.FontSize = 20;
    grid on;

    figure;
    scatter(E_min, test_score, 'filled', 'b');
    xlabel('Fuzzy Entropy', 'FontSize', 25);
    ylabel('Sleep Efficiency', 'FontSize', 25);
    title(['Scatter Plot at Min Correlation Time Scale ', num2str((min_idx + num_2) * 5)], 'FontSize', 25);
    ax = gca;
    ax.FontSize = 20;
    grid on;

    disp(['Max Correlation Time Scale: ', num2str((max_idx + num_2) * 5)]);
    disp(['Min Correlation Time Scale: ', num2str((min_idx + num_2) * 5)]);
end
% % %%%%%%%%%%%%%%%%%%%%%
% % scale sample
%     figure;
%     hold on;
%     %time_scales = (scale_bottom+num_2) * 5:5: (scale+num_2) * 5;
%     significant = p_values < 0.05;
%     non_significant = ~significant;
% 
%     scatter((scale_bottom+num_2):(scale+num_2),correlation_coefficients(significant), 'filled', 'b'); %r
%     scatter((scale_bottom+num_2):(scale+num_2),correlation_coefficients(non_significant), 'filled', 'b');
% 
%     % 境界線の追加
%     yline(min(correlation_coefficients(significant)), '--r', 'p = 0.05', 'LabelHorizontalAlignment', 'left','LineWidth',5);
% 
%     xlabel('Time Scale [s]');
%     ylabel('Spearman Correlation Coefficient');
%     title('Correlation Coefficients over Different Time Scales');
%     %legend('p < 0.05', 'p \geq 0.05');
%     grid on;
%     hold off;
% 
% end
