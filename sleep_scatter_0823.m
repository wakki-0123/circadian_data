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

function [p_values, correlation_coefficients, t_values] = sleep_scatter_0823(e)
    % 被験者のテストスコア
    % test_score = [0.886324022, 0.866727644, 0.911264649, 0.877716281, 0.885269186, ...
    %               0.855124121, 0.871734582, 0.859731051, 0.900972038, 0.854969546, ...
    %               0.8519555, 0.894203152,0.905361808];

    test_score = [0.876847548,0.867301668,0.888202161,0.87524095,0.85704005,0.911264649,0.855918089,0.901118598,0.859571395,0.886905742,0.8519555,0.916248705];
 
    entropy = e; % 12x10000 行列
    num = 12; % 被験者の数
    scale = 200; % スケール数を変更

    % 相関係数、p値、t値の初期化
    correlation_coefficients = zeros(1, scale);
    p_values = zeros(1, scale);
    t_values = zeros(1, scale);

    % 各時間スケールのデータを処理
    for i = 1:scale
        % 各時間スケールごとにエントロピー値を取り出す
        E_person = entropy(1:num, i); % numx1 ベクトル
        
        % Spearmanの順位相関係数を計算
        [correlation_coefficients(i), p_values(i)] = corr(E_person, test_score(1:num)', 'Type', 'Spearman');
        
        % t値の計算
        n = num; % サンプル数
        t_values(i) = correlation_coefficients(i) * sqrt((n - 2) / (1 - correlation_coefficients(i)^2));
    end

    % 相関係数の図の作成
    figure;
    scatter(1:scale, correlation_coefficients, 'filled');
    xlabel('Time Scale');
    ylabel('Spearman Correlation Coefficient');
    title('Correlation Coefficients over Different Time Scales');
    grid on;
end
