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


errp_original_plus=mean(e_all_array_0411)+std(e_all_array_0411);
errp_original_minus=mean(e_all_array_0411)-std(e_all_array_0411);

errp_IAAFT_plus = mean(e_all_array_IAAFT_0411) + std(e_all_array_IAAFT_0411);
errp_IAAFT_minus = mean(e_all_array_IAAFT_0411) - std(e_all_array_IAAFT_0411);




% グラフの表示
figure;

%plot(time, e1, 'r')
% errorbar(time,mean(e_all_array_0411),std(e_all_array_0411),'--r');
% hold on
% errorbar(time, mean(e_all_array_IAAFT_0411), std(e_all_array_IAAFT_0411), '--b');
% legend('ORG', 'IAAFT', 'Location', 'southeast');
loglog(time,mean(e_all_array_0411),'r','LineWidth',5)
hold on
loglog(time,errp_original_plus,'--r','LineWidth',3)
loglog(time,errp_original_minus,'--r','LineWidth',3)

loglog(time,mean(e_all_array_IAAFT_0411),'b','LineWidth',5)
loglog(time,errp_IAAFT_plus,'--b','LineWidth',3)
loglog(time,errp_IAAFT_minus,'--b','LineWidth',3)
lgd = legend('ORG', 'ORG ERRORBAR','IAAFT ERRORBAR','IAAFT', 'Location', 'northeast');
lgd.FontSize = 40;
%set(gca, 'XScale', 'log');
hold off
ax = gca;
ax.FontSize = 40;
grid on;
ylim([0.1 3]);
%title('Heart Rate Multiscale Fuzzy Entropy');
xlabel('Time Scale [sec]');
ylabel('Fuzzy Entropy');
