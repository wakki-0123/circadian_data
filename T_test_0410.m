function [p,t] = T_test_0410(data1, data2)
    

    % 結果を格納する変数を初期化
    % num_iterations = length(data1);  % data1とdata2は同じ長さと仮定
    % p = NaN(1, num_iterations);
    % t = NaN(1, num_iterations);
    disp(data1(1:11))
    alpha = 0.05;

    % 各データ点のペアに対してt検定を実行
    for i = 1:10000

        [h, p(i), ci, stats] = ttest(data1(:,i), data2(:,i));
        t(i) = stats.tstat;
    end
    % p値を用いて帰無仮説の棄却/採択を判断

% p_corrected = mafdr(p, 'BHFDR', true);
% fprintf('補正されたp値: %.4f\n', p_corrected);




    

figure

plot(t)







