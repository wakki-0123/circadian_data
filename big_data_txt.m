function big_data_txt()
% テキストファイルを読み込む
filename = 'heartrate_log_20240727-20240816_000000000000.txt'; % ここにあなたのファイル名を指定
data = readtable(filename);

% user_id が「2」のデータをフィルタリング
%filteredData = data(strcmp(data.user_id, '002'), :);
filteredData = data(data.user_id == 3, :);


% 結果を表示または保存
disp(filteredData);
writetable(filteredData, 'kawakami_heartrate2.xlsx'); % フィルタされたデータを別のファイルに保存
