function [e_all, e_IAAFT_all] = MFE_circadian_multiple(data_cell, c, maxiter, m, factor, mf, rn, local, tau)
% MFE_circadian_multiple: Circadian Rhythm用の複数データを処理する関数

% 0405に実行したプログラム
% MFE_circadian_multiple(data_cell,10,50,2,10000,'Exponential',[0.2 2],0,1);

% データの数
num_data = numel(data_cell);

% 出力用の変数を初期化
e_all = cell(1, num_data);
e_IAAFT_all = cell(1, num_data);
q = 0;

for data_index = 1:2
    % 各データを取得
    q = q + 1;
    data = data_cell{data_index};
    data_l = length(data);
    
    % マルチスケールファジーエントロピーとIAAFTを計算
    [e, e_IAAFT] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau, data_l);
    
    % 結果を保存
    e_all{data_index} = e;
    e_IAAFT_all{data_index} = e_IAAFT;
    disp("現在の稼働状況10段階 (全体でどのくらい進んでいるか)")
    disp(q)
end

end

function [e1, e_IAAFT] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau, data_l)
% MFE_circadian_single: 単一のデータに対してマルチスケールファジーエントロピーとIAAFTを計算する関数

% マルチスケールファジーエントロピーの計算
e1 = fuzzymsentropy(data, m, mf, rn, local, tau, factor);

% IAAFTを実行し、その結果のマルチスケールファジーエントロピーを計算
e2 = zeros(factor, c);
[s, ~] = IAAFT(data, c, maxiter);

% サロゲートデータのマルチスケールファジーエントロピーを計算
for i = 1:c
    e2(:, i) = fuzzymsentropy(s(:, i), m, mf, rn, local, tau, factor);
    disp("現在のサロゲートデータに関するファジーエントロピーの様子")
    disp(i)
end

% 転置行列
e2 = e2';
e_IAAFT = mean(e2);

% グラフの表示
plot_MFE_graph(e1, e2, data_l);

end

function plot_MFE_graph(e1, e2, data_l)
% マルチスケールファジーエントロピーのグラフをプロットする関数
time_length = data_l * 5; % 全部の区間の秒数
factor = size(e2, 2);

% 時間スケールを計算
time_s = zeros(1, factor);
time = zeros(1, factor);
for i = 1:factor
    time_s(i) = data_l / i; % 合計サンプルの個数
    time(i) = time_length / time_s(i); % タイムスケール
end

% グラフの表示
figure;
plot(time, e1, 'r')
hold on
errorbar(time, mean(e2), std(e2), 'b');
legend('ORG', 'IAAFT', 'Location', 'southeast');
set(gca, 'XScale', 'log');
hold off
title('Heart Rate Multiscale Fuzzy Entropy');
xlabel('Time Scale');
ylabel('Fuzzy Entropy');

end
%%%%%%%%%%%%%%%%%%%%%%%%%
% 0627 add
% function e = fuzzymsentropy(input, m, mf, rn, local, tau, factor)
%     y = input;
%     y = (y - mean(y)) / std(y);
% 
%     e = zeros(1, factor);
%     for i = 1:factor
%         if i == 1
%             s = y;
%         else
%             s = coarsegraining(y, i);
%         end
%         e(i) = FuzEn_MFs(s, m, mf, rn, local, tau);
%     end
%     e = e';
% end

