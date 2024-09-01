% function sleep_soukan(e_all_array_0411)
%     % 被験者のテストスコア
%     test_score = [0.876847548,0.867301668,0.888202161,0.87524095,0.85704005,0.911264649,0.855918089,0.901118598,0.859571395,0.886905742,0.8519555, 0.894203152];
% 
%     % データのサイズ取得
%     num_subjects = 12;
%     num_time_scales = 200;
% 
% 
%     % 各被験者ごとのプロットを作成
%     figure;
%     for i = 1:num_subjects
%         % i番目の被験者のデータとテストスコアの関係をプロット
%         subplot(4, 3, i);  % 4行3列のサブプロットに配置
%         plot(1:num_time_scales, e_all_array_0411(i, 1:num_time_scales), '-o');
%         title(['Subject ', num2str(i)]);
%         xlabel('Time Scales');
%         ylabel('Data Value');
%         ylim([0.1 1]);
%         grid on;
% 
%         % 縦軸にテストスコアを表示
%         yyaxis right
%         plot(1:num_time_scales, repmat(test_score(i), 1, num_time_scales), '--r');
%         ylabel('Test Score');
%     end
% 
%     % グラフ全体のタイトル
%     sgtitle('Data vs. Test Score across Time Scales');
% end
function sleep_soukan(e_all_array_0411)
    % 被験者のテストスコア
    test_score = [0.876847548,0.867301668,0.888202161,0.87524095,0.85704005,0.911264649,0.855918089,0.901118598,0.859571395,0.886905742,0.890027395,0.916248705]';

    % データのサイズ取得
    num_subjects = 13;
    num_time_scales = 200;

    

% 相関係数の計算
correlations = zeros(1, 200); % 相関係数を格納するためのベクトル

for i = 1:num_time_scales
    correlations(i) = corr(test_score(1:10), e_all_array_0411(1:10, i), 'Type', 'Spearman'); % スペアマンの相関係数を計算
end

% 相関係数を降順にソートし、上位10個を取得
    [sorted_corr, sorted_idx] = sort(correlations, 'descend');

    % 結果の表示 - 上位10個
    fprintf('上位10個の最も強い相関を持つデータセット:\n');
    for j = 1:10
        fprintf('第 %d 位: 列 %d, 相関係数 = %f\n', j, sorted_idx(j), sorted_corr(j));
    end
end
