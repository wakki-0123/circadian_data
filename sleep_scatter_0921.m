function [p_values, correlation_coefficients, t_values, peak_info] = sleep_scatter_0921(e)
    % 被験者のテストスコア
    test_score = [0.876847548,0.867301668,0.888202161,0.87524095,0.85704005,0.911264649,0.855918089,0.901118598,0.859571395,0.886905742,0.891167085,0.8519555,0.890027395,0.890409196,0.916248705];

    entropy = e; % 15x10000 行列
    num = 15; % 被験者の数
    scale = 2000; % スケール数を変更
    scale_bottom = 1; % スケールの底

    % 相関係数、p値、t値の初期化
    correlation_coefficients = zeros(1, scale-scale_bottom+1);
    p_values = zeros(1, scale-scale_bottom+1);
    t_values = zeros(1, scale-scale_bottom+1);
    num_2 = 7;
    
    % 各時間スケールのデータを処理
    for i = scale_bottom:1:scale
        % 各時間スケールごとにエントロピー値を取り出す
        E_person = entropy(1:num, i); % numx1 ベクトル
        
        % Spearmanの順位相関係数を計算
        [correlation_coefficients(i - scale_bottom + 1), p_values(i - scale_bottom + 1)] = corr(E_person, test_score(1:num)', 'Type', 'Spearman');
        
        % t値の計算
        n = num; % サンプル数
        t_values(i - scale_bottom + 1) = correlation_coefficients(i - scale_bottom + 1) * sqrt((n - 2) / (1 - correlation_coefficients(i - scale_bottom + 1)^2));
    end

    % スムージング処理（例：Savitzky-Golay フィルタ）
    smooth_coefficients = sgolayfilt(correlation_coefficients, 3, 11);

    % ピークの検出
    [peak_value, peak_index] = max(smooth_coefficients);
    peak_scale = peak_index + scale_bottom - 1;
    
    % ピーク情報の保存
    peak_info.peak_value = peak_value;
    peak_info.peak_scale = peak_scale;

    % ホルモン動態モデルの例
    time = linspace(scale_bottom, scale, scale-scale_bottom+1); % 相関係数データと同じスケール範囲に合わせる
    cortisol = exp(-((time - (scale/2)).^2) / 40000); % コルチゾールの逆U字型分泌パターン

    % 自律神経系モデルの例
    time_hrv = 0:0.1:24; % 24時間のスケール
    hrv = cos(time_hrv / 2) .* exp(-((time_hrv - 12).^2) / 4); % HRVの逆U字型パターン

    % 相関係数の図の作成
    figure;
    plot(time, correlation_coefficients, 'o', 'DisplayName', 'Correlation Coefficients');
    hold on;
    plot(time, smooth_coefficients, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed Data');
    plot(peak_scale, peak_value, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b', 'DisplayName', 'Peak');
    plot(time, cortisol, 'g--', 'LineWidth', 2, 'DisplayName', 'Cortisol Pattern');
    plot(linspace(scale_bottom, scale, length(time_hrv)), hrv, 'm-', 'LineWidth', 2, 'DisplayName', 'HRV Pattern');
    xlabel('Factor Scale');
    ylabel('Spearman Correlation Coefficient');
    title('Correlation Coefficients over Different Time Scales and Physiological Patterns');
    legend;
    grid on;
end
