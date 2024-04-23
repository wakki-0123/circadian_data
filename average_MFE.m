function [e1,e2] = average_MFE(e_all_array_0411,e_all_array_IAAFT_0411)
data_l = 10000;
time_length = data_l * 5; % 全部の区間の秒数
factor = 10000;

% 時間スケールを計算
time_s = zeros(1, factor);
time = zeros(1, factor);
for i = 1:factor
    time_s(i) = data_l / i; % 合計サンプルの個数
    time(i) = time_length / time_s(i); % タイムスケール
end
e1 = mean(e_all_array_0411);



% グラフの表示
figure;
plot(time, e1, 'r')
hold on
errorbar(time, mean(e_all_array_IAAFT_0411), std(e_all_array_IAAFT_0411), 'b');
legend('ORG', 'IAAFT', 'Location', 'southeast');
set(gca, 'XScale', 'log');
hold off
title('Heart Rate Multiscale Fuzzy Entropy');
xlabel('Time Scale');
ylabel('Fuzzy Entropy');
