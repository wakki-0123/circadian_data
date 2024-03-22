% % 対標本t検定を実行
% [h, p, ci, stats]=ttest(person_ans_heartrate_AllDay, person_IAAFT_MSEaverage_AllDay);
% 
% % 結果の表示
% fprintf('t統計量: %.4f\n', stats.tstat);
% fprintf('p値: %.4f\n', p);
% 
% % 有意水準を設定
% alpha = 0.05;
% 
% % p値を用いて帰無仮説の棄却/採択を判断
% if p < alpha
%     fprintf('帰無仮説を棄却: サンプル間に統計的に有意な違いがあります\n');
% else
%     fprintf('帰無仮説を採択: サンプル間に統計的に有意な違いがありません\n');
% end



% 有意水準を設定
alpha = 0.05;

% ループを使用してデータごとに対標本t検定を実行
for i = 1:12000
%     disp(['MSE: ',num2str(person_ans_heartrate_AllDay_1(i))]);
%     disp(['IAAFTafMSE: ',num2str(person_IAAFT_MSEaverage_AllDay_1(i))]);    
% 
%     sample1 = person_ans_heartrate_AllDay_1(i);  % サンプル1から1つのデータを選択
%     sample2 = person_IAAFT_MSEaverage_AllDay_1(i);  % サンプル2から対応するデータを選択

    x=person_ans_heartrate_AllDay_1(:,1);
    y=person_IAAFT_MSEaverage_AllDay_1(:,1);

    disp(x)

    [h, p(i), ci, stats] = ttest2(x, y);

    % 対標本t検定を実行
%     [h, p, ci, stats] = ttest2(sample1, sample2);
    

    %結果の表示

   % fprintf('データセット %d に対するt統計量: %.4f\n', i,stats.tstat);
    %fprintf('データセット %d に対するt統計量: %.4f, p値: %.4f\n', i, stats.tstat, p);
% 
%     % p値を用いて帰無仮説の棄却/採択を判断
%     if p < alpha
%         fprintf('データセット %d: 帰無仮説を棄却, 統計的に有意な違いがあります\n', i);
%     else
%         fprintf('データセット %d: 帰無仮説を採択, 統計的に有意な違いがありません\n', i);
%     end
end