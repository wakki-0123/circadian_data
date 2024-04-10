function threshold = calculate_FDR_threshold(p_values, alpha)
    % p_values: データ点に対するp値の配列
    % alpha: FDRの許容誤り率（例えば0.05）

    % FDR補正を適用
    [p_sorted, idx] = sort(p_values);
    m = length(p_values);  % データ点の総数
    adjusted_alpha = alpha * (1:m) / m;

    % FDR補正後のしきい値の計算
    significant_idx = find(p_sorted <= adjusted_alpha);
    if isempty(significant_idx)
        threshold = 0;  % 有意なデータ点が見つからない場合
    else
        max_significant_idx = max(significant_idx);
        threshold = p_sorted(max_significant_idx);
    end
end