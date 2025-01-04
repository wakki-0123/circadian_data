function sleep_soukan(e_all_array_0411)
    % 被験者のテストスコア
    test_scores = [0.876847548, 0.867301668, 0.888202161, 0.87524095, 0.85704005, ...
                   0.911264649, 0.855918089, 0.901118598, 0.859571395, 0.886905742, ...
                   0.891167085, 0.8519555, 0.890027395, 0.890409196, 0.916248705];

    % データのサイズ取得
    num_subjects = length(test_scores);
    [num_rows] = size(e_all_array_0411);
    num_time_scales = 400;

    % データの行数チェック
    if num_rows < num_subjects
        error('e_all_array_0411 の行数は %d 以上でなければなりません。', num_subjects);
    end

    % 相関係数の計算
    correlations = zeros(1, num_time_scales); % 相関係数を格納するためのベクトル

    for i = 1:num_time_scales
        correlations(i) = corr( e_all_array_0411(1:num_subjects, i),test_scores', 'Type', 'Spearman'); % スペアマンの相関係数を計算
    end

    % 相関係数を降順にソートし、上位10個を取得
    [sorted_corr, sorted_idx] = sort(correlations, 'descend');

    % 結果の表示 - 上位10個
    fprintf('上位10個の最も強い相関を持つデータセット:\n');
    for j = 1:10
        fprintf('第 %d 位: 列 %d, 相関係数 = %f\n', j, sorted_idx(j), sorted_corr(j));
    end
end

