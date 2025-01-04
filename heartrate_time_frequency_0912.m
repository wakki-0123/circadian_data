function heartrate_time_frequency_0912()
    % ファイルパスを指定
    file = ['C:/Users/iw200/Desktop/NEW_heartrate_kari/wakita_heartrate.xlsx'];
    
    % Excelファイルからデータをテーブル形式で読み込む
    T = readtable(file);

    % デバッグ用に読み込んだテーブルデータを表示
    %disp('Excelファイルから読み込んだデータ:');
    %disp(T);

    % 日時データを抽出
    datetimeData = T.datetime; 

    % datetimeDataをMATLABのdatetime形式に変換
    y1 = datetime(datetimeData, 'InputFormat', 'yyyy/MM/dd HH:mm');

    % NaT値があるかどうかをチェックし、メッセージを表示
    if any(isnat(y1))
        disp('警告: 一部の日時値が変換できずNaTになっています。');
    end

    % 日付ごとにデータをグループ化
    [uniqueDates, ~, idx] = unique(dateshift(y1, 'start', 'day'));
    numDays = length(uniqueDates);

    % 日ごとのサンプリング周波数を計算
    samplingFrequencyPerDay = zeros(numDays, 1);
    for i = 1:numDays
        % 当日の日時データを取得
        dayData = y1(idx == i);
        
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
