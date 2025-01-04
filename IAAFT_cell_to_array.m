function averaged_data = IAAFT_cell_to_array(data)

for i = 1:3
    averaged_data(i, :) = mean(data{i}, 1); % 行ごとに平均を計算
end
