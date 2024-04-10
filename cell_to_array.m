function [values,means1] = cell_to_array(e_all,e_IAAFT_all)
% 要素数を取得
num_elements = numel(e_all{1});
num_elements1 = numel(e_IAAFT_all{1});
 % 平均値を格納する配列を作成
means1 = zeros(1, num_elements1); 

for i = 1:num_elements
    % 各要素番号ごとの値を取得
    values = cellfun(@(x) x(i), e_all);
    values1 = cellfun(@(x) x(i), e_IAAFT_all);

    % 平均を計算して格納
    %means(i) = mean(values);
    means1(i) = mean(values1);

end

% 結果を表示
disp('平均値:');
disp(values);
disp(means1);