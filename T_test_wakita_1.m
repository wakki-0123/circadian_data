function p1 = T_test_wakita_1(person_ans_heartrate_AllDay, person_IAAFT_MSEaverage_AllDay)
% 対標本t検定を実行


[h, p, ci, stats]=ttest(person_ans_heartrate_AllDay, person_IAAFT_MSEaverage_AllDay);

% 結果の表示
fprintf('t統計量: %.4f\n', stats.tstat);
fprintf('p値: %.4f\n', p);
p_corrected = mafdr(p, 'BHFDR', true);
fprintf('補正されたp値: %.4f\n', p_corrected);



% 有意水準を設定
alpha = 0.05;

% p値を用いて帰無仮説の棄却/採択を判断
if p_corrected < alpha
    fprintf('帰無仮説を棄却: サンプル間に統計的に有意な違いがあります\n');
else
    fprintf('帰無仮説を採択: サンプル間に統計的に有意な違いがありません\n');
end
p1 = p_corrected;