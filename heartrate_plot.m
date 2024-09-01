function y1 = heartrate_plot()

% 心拍変動を表示するための関数
file = 'C:/Users/iw200/Desktop/NEW_heartrate_kari/asato_heartrate.xlsx';
%データの読み込み
% file: 任意の被験者のファイル名

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excelファイルから数値データを読み込む
[numData, ~, ~] = xlsread(file);

% データの確認
%disp(numData(:,3));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = numData(:,3);
figure;
plot(y)

legend('ORG','Location','southeast');
title('original heartrate data');
xlabel('sample');
ylabel('heartrate');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%前処理
l = length(y);
t_original = 1:1:l;
% 新しいサンプリング周波数を指定
new_sampling_frequency = 0.2;  % 新しいサンプリング周波数 (Hz)

% 新しい時間ベクトルを自動的に計算
new_time_vector = 1:1/new_sampling_frequency:max(t_original);

% resample関数を使用してデータをリサンプリング
y_resampled = interp1(t_original, y, new_time_vector);

y1 = y_resampled;
y1 = zscore(y1); % Zscore化

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 前処理後の時系列データのplot
figure;
plot(y1)

legend('ORG','Location','southeast');
title('resample heartrate data');
xlabel('time[s]');
ylabel('heartrate');

disp(length(y1))
