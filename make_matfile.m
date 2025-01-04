function make_matfile(data_cell)
% 大きなデータを生成


% matfile でデータを保存
filename = 'heart_data_0906.mat';
matObj = matfile(filename, 'Writable', true);
matObj.data = data_cell;
