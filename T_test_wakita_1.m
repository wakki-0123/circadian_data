function p1 = T_test_wakita_1(person_ans_heartrate_AllDay, person_IAAFT_MSEaverage_AllDay)
% 対標本t検定を実行


[h, p, ci, stats]=ttest(person_ans_heartrate_AllDay, person_IAAFT_MSEaverage_AllDay);

% 結果の表示
fprintf('t統計量: %.4f\n', stats.tstat);
fprintf('p値: %.4f\n', p);
p_corrected = mafdr(p, 'BHFDR', true);
fprintf('補正されたp値: %.4f\n', p_corrected);



%有意水準を設定
alpha = 0.05;

% p値を用いて帰無仮説の棄却/採択を判断
if p < alpha
    fprintf('帰無仮説を棄却: サンプル間に統計的に有意な違いがあります\n');
else
    fprintf('帰無仮説を採択: サンプル間に統計的に有意な違いがありません\n');
end
% ループを使用してデータごとに対標本t検定を実行
for i = 1:10000 % 前は12000
%     disp(['MSE: ',num2str(person_ans_heartrate_AllDay_1(i))]);
%     disp(['IAAFTafMSE: ',num2str(person_IAAFT_MSEaverage_AllDay_1(i))]);    
% 
%     sample1 = person_ans_heartrate_AllDay_1(i);  % サンプル1から1つのデータを選択
%     sample2 = person_IAAFT_MSEaverage_AllDay_1(i);  % サンプル2から対応するデータを選択

    x=person_ans_heartrate_AllDay(:,1);
    y=person_IAAFT_MSEaverage_AllDay(:,1);


    [h, p, ci, stats]=ttest(x,y);

  p1(i) = p;
  t(i) = stats.tstat;
end

figure
plot(p1)
hold on
plot(t)



