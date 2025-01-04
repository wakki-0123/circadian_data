function sleep_peak(e_all_array_0411)
    % 被験者のテストスコア 13 , 12人
     test_score = [0.876847548,0.867301668,0.888202161,0.87524095,0.85704005,0.911264649,0.855918089,0.901118598,0.859571395,0.886905742,0.891167085,0.8519555,0.890027395,0.890409196,0.916248705]';
 
    %test_score = [0.876847548,0.867301668,0.888202161,0.87524095,0.85704005,0.911264649,0.855918089,0.901118598,0.859571395,0.886905742,0.8519555,0.916248705]';

    % データのサイズ取得
    num_subjects = 15;
    num_time_scales = 200;

    % 相関係数の計算
    correlations = zeros(1, num_time_scales); % 相関係数を格納するためのベクトル

    for i = 40:50
        figure; % 新しい図を作成

        % スペアマンの相関係数を計算
        correlations(i) = corr(test_score, e_all_array_0411(1:num_subjects, i), 'Type', 'Spearman');
        
        % エントロピー値を横軸に、テストスコアを縦軸にプロット
        plot(e_all_array_0411(:, i), test_score, 'o-', 'DisplayName', ['Time Scale ' num2str(i)]);
        
        % プロットの設定
        xlabel('エントロピー値');
        ylabel('テストスコア');
        title(['Time Scale ' num2str(i) ': エントロピー値とテストスコアのプロット']);
        legend('show'); % 凡例を表示
    end

    % 相関係数を降順にソートし、上位10個を取得
    [sorted_corr, sorted_idx] = sort(correlations, 'descend');

    % 結果の表示 - 上位10個
    fprintf('上位10個の最も強い相関を持つデータセット:\n');
    for j = 1:10
        fprintf('第 %d 位: 列 %d, 相関係数 = %f\n', j, sorted_idx(j), sorted_corr(j));
    end
end
