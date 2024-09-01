function kari_kaiseki(correlation_coefficients)

[pks, locs] = findpeaks(correlation_coefficients, 'MinPeakProminence', 0.01);
peakScale = locs(pks == max(pks)); % 最も高いピークを特定
range = 10; % ピークの前後10スケール
startIdx = max(1, peakScale - range);
endIdx = min(length(correlation_coefficients), peakScale + range);

% ピーク周辺のデータを抽出
detailed_correlation = correlation_coefficients(startIdx:endIdx);

% ピーク周辺のプロット
figure;
plot(startIdx:endIdx, detailed_correlation, 'o-', 'LineWidth', 2);
xlabel('Time Scale');
ylabel('Spearman Correlation Coefficient');
title('Detailed Analysis around the Peak');
grid on;
fit_curve = polyfit(startIdx:endIdx, detailed_correlation, 2);
fitted_values = polyval(fit_curve, startIdx:endIdx);

% フィッティング結果のプロット
figure;
plot(startIdx:endIdx, detailed_correlation, 'o');
hold on;
plot(startIdx:endIdx, fitted_values, '-r');
xlabel('Time Scale');
ylabel('Correlation Coefficient');
title('Fitting around the Peak');
grid on;
hold off;

numBootstrap = 1000;
bootstrap_means = bootstrp(numBootstrap, @mean, detailed_correlation);

% ブートストラップの結果をヒストグラムで表示
figure;
histogram(bootstrap_means);
xlabel('Mean Correlation Coefficient');
ylabel('Frequency');
title('Bootstrap Analysis around the Peak');
grid on;

% 相関係数の1次微分（変化率）
first_derivative = diff(correlation_coefficients);

% 相関係数の2次微分（曲率）
second_derivative = diff(first_derivative);

% 微分のプロット
figure;
subplot(2,1,1);
plot(first_derivative, 'LineWidth', 2);
xlabel('Time Scale');
ylabel('First Derivative');
title('First Derivative of Correlation Coefficients');

subplot(2,1,2);
plot(second_derivative, 'LineWidth', 2);
xlabel('Time Scale');
ylabel('Second Derivative');
title('Second Derivative of Correlation Coefficients');
grid on;
% ウェーブレット変換の実行
[c, l] = wavedec(correlation_coefficients, 4, 'db1'); % db1はDaubechiesのウェーブレット

% ピーク周辺のスケール範囲で逆変換
peak_coefficients = wrcoef('a', c, l, 'db1', 4);

% 結果のプロット
figure;
plot(peak_coefficients);
xlabel('Time Scale');
ylabel('Wavelet Coefficient');
title('Wavelet Coefficients around the Peak');
grid on;

% ピーク周辺のZスコア計算
mean_val = mean(correlation_coefficients);
std_val = std(correlation_coefficients);
z_scores = (correlation_coefficients - mean_val) / std_val;

% 異常値のプロット
figure;
plot(z_scores, 'o-');
xlabel('Time Scale');
ylabel('Z-Score');
title('Z-Score Analysis around the Peak');
grid on;


