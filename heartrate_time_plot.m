function y1 = heartrate_time_plot()
% 0926現在、心拍数データのリサンプリングはこのプログラムを使用している

% ファイルパスを指定

file = 'C:/Users/iw200/Desktop/NEW_heartrate_kari/sasaki_heartrate.xlsx';
%file = 'ebato_heartrate.xlsx';

% Excelファイルからデータをテーブル形式で読み込む
T = readtable(file);

% 日時データと心拍数データを抽出
datetimeData = T.datetime;
heartRateData = T.heartrate;

% datetimeDataをMATLABのdatetime形式に変換
time = datetime(datetimeData, 'InputFormat', 'yyyy/MM/dd HH:mm');

% NaT値があるかどうかをチェックし、メッセージを表示
if any(isnat(time))
    disp('警告: 一部の日時値が変換できずNaTになっています。');
end

% 元のサンプリング周波数を計算
sampling_intervals = seconds(diff(time));
original_sampling_frequency = 1 / mean(sampling_intervals);
disp(original_sampling_frequency)

% 新しいサンプリング周波数を指定
new_sampling_frequency = 0.2;  % 新しいサンプリング周波数 (Hz)

% リサンプリング係数の計算
% サンプリング周波数の比率を整数に近似
resample_ratio = original_sampling_frequency / new_sampling_frequency;

% 最も近い整数の分子と分母を求める
[p, q] = rat(resample_ratio);
disp(p);
disp(q);

% データのリサンプリング
% resample関数を使用してリサンプリング
y_resampled = resample(heartRateData, q, p);

% 新しいサンプル数の計算
new_num_samples = length(y_resampled);

% Z-score化
y1 = zscore(y_resampled);

% 新しい時間ベクトルの生成
% 新しいサンプル数を使用して新しい時間ベクトルを生成
new_time_vector = linspace(time(1), time(end), new_num_samples);

% リサンプリング後のサンプリング間隔とサンプリング周波数を計算
resampled_time_intervals = diff(new_time_vector);
resampled_time_intervals_seconds = seconds(resampled_time_intervals);
resampled_sampling_frequency = 1 / mean(resampled_time_intervals_seconds);

% リサンプリング後の時系列データのプロット
figure;
plot(new_time_vector, y1)
legend('Resampled and Z-scored','Location','southeast');
title('Resampled heartrate data');
xlabel('Time');
ylabel('Heartrate');

% 結果を表示
disp(length(y1))
disp('Original Sampling Frequency:');
disp(original_sampling_frequency)
disp('Resampled Sampling Frequency:');
disp(resampled_sampling_frequency)

% 日付ごとにデータをグループ化
[uniqueDates, ~, idx] = unique(dateshift(time, 'start', 'day'));
numDays = length(uniqueDates);

% 日ごとのサンプリング周波数を計算
samplingFrequencyPerDay = zeros(numDays, 1);
for i = 1:numDays
    % 当日の日時データを取得
    dayData = time(idx == i);
    
    % 当日の時間差を計算
    timeDifferences = diff(dayData);
    
    % 当日のサンプリング周波数を計算 (1秒の時間差で)
    if ~isempty(timeDifferences)
        samplingFrequencyPerDay(i) = mean(1 ./ seconds(timeDifferences));
    else
        samplingFrequencyPerDay(i) = NaN;
    end
end

% 結果を表示
results = table(uniqueDates, samplingFrequencyPerDay, 'VariableNames', {'Date', 'SamplingFrequency'});
disp(results);

end

