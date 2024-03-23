function y1 = heartrate_plot(file)
%file = 'fitbit_heartrate_修論用_安田.xlsx';
%データの読み込み

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%file = 'before_pupil_diameter.csv'; %任意のファイルのパス
%file1 = 'after_pupil_diameter.csv'
% CSVファイルを読み取る

% CSVファイルを読み取る
% Excelファイルから数値データを読み込む
[numData, ~, ~] = xlsread(file);

% データの確認
%disp(numData(:,3));





%前処理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = numData(:,3);
%y = medfilt1(data,thredshold,'omitnan','truncate');%メジアン補間(NaNの処理)
%y = fillmissing(y,'linear');%線形補間
%y = filloutliers(y,'linear','movmedian',thredshold);



figure;
plot(y)

legend('ORG','Location','southeast');
title('original heartrate data');
xlabel('sample');
ylabel('heartrate');

l = length(y);
t_original = 1:1:l;
% 新しいサンプリング周波数を指定
new_sampling_frequency = 0.2;  % 新しいサンプリング周波数 (Hz)

% 新しい時間ベクトルを自動的に計算
new_time_vector = 1:1/new_sampling_frequency:max(t_original);

% resample関数を使用してデータをリサンプリング
y_resampled = interp1(t_original, y, new_time_vector);



y1 = y_resampled;
y1 = zscore(y1);


%y = data; %元のデータを表示したいときはコメントを外す(欠損値だらけ)
% 前処理後の時系列データのplot
figure;
plot(y1)

legend('ORG','Location','southeast');
title('resample heartrate data');
xlabel('time[s]');
ylabel('heartrate');
