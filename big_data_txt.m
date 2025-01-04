function big_data_txt()
% テキストファイルを読み込む
filename = 'heartrate_log_20241111-20241201_000000000000.csv'; % ここにあなたのファイル名を指定
data = readtable(filename);

% user_id が「2」のデータをフィルタリング
%filteredData = data(strcmp(data.user_id, '002'), :);
filteredData = data(data.user_id == 11, :);


% 結果を表示または保存
disp(filteredData);
writetable(filteredData, 'sasaki_heartrate.xlsx'); % フィルタされたデータを別のファイルに保存
