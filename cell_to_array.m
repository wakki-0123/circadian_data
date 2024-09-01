function [values] = cell_to_array(e_all)
    % 各セルの中にある配列の長さを取得
    num_elements = numel(e_all{1});
    
    % 行列の初期化 (行数はセル配列の数、列数は各配列の長さ)
    num_cells = numel(e_all);
    values = zeros(num_cells, num_elements);
    
    % 各要素を行列に格納
    for i = 1:num_elements
        % 各セル配列の i 番目の要素を行列の i 番目の列に格納
        values(:, i) = cellfun(@(x) x(i), e_all);
    end
end
