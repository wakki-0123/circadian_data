function [p,t] = T_test_0410(data1, data2)

% 1210
% 結果を格納する変数を初期化
% num_iterations = length(data1);  % data1とdata2は同じ長さと仮定
% p = NaN(1, num_iterations);
% t = NaN(1, num_iterations);

alpha = 0.05;
%data_l = 9993;
data_l = size(data1, 2);
% 各データ点のペアに対してt検定を実行
for i = 1:data_l

    [h, p(i), ci, stats] = ttest(data1(:,i), data2(:,i));
    t(i) = stats.tstat;

    % p値を用いて帰無仮説の棄却/採択を判断
    if p < alpha
        fprintf('帰無仮説を棄却: サンプル間に統計的に有意な違いがあります\n');
    else
        fprintf('帰無仮説を採択: サンプル間に統計的に有意な違いがありません\n');
    end
end


figure

plot(t)







