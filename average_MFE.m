function average_MFE(e_all_array_0411,e_all_array_IAAFT_0411)
% ただファジイエントロピーをplotするだけ 1210
data_l = 9993;
data_gyou = size(e_all_array_0411,1);
num = 7; % とばしたタイムスケール
time_length = (data_l + num) * 5; % 全部の区間の秒数
factor = 9993;

% 時間スケールを計算
time_s = zeros(1, factor);
time = zeros(1, factor);
for i = 1:factor
    time_s(i) = data_l / i; % 合計サンプルの個数
    time(i) = time_length / time_s(i); % タイムスケール
end

e1 = mean(e_all_array_0411,1);
e2 = mean(e_all_array_IAAFT_0411,1);

%%%%%%%%%%%
% STE version
errp_original_plus=e1+(std(e_all_array_0411,1)/sqrt(data_gyou));
errp_original_minus=e1-(std(e_all_array_0411,1)/sqrt(data_gyou));

errp_IAAFT_plus = e2 + (std(e_all_array_IAAFT_0411,1)/sqrt(data_gyou));
errp_IAAFT_minus = e2 - (std(e_all_array_IAAFT_0411,1)/sqrt(data_gyou));


%%%%%
% STD version
% errp_original_plus=e1+std(e_all_array_0411,1);
% errp_original_minus=e1-std(e_all_array_0411,1);
% 
% errp_IAAFT_plus = e2 + std(e_all_array_IAAFT_0411,1);
% errp_IAAFT_minus = e2 - std(e_all_array_IAAFT_0411,1);

% グラフの表示
figure;

%plot(time, e1, 'r')
% errorbar(time,mean(e_all_array_0411),std(e_all_array_0411),'--r');
% hold on
% errorbar(time, mean(e_all_array_IAAFT_0411), std(e_all_array_IAAFT_0411), '--b');
% legend('ORG', 'IAAFT', 'Location', 'southeast');
loglog(time,e1,'r','LineWidth',2)
hold on
loglog(time,errp_original_plus,'--r','LineWidth',3)
loglog(time,errp_original_minus,'--r','LineWidth',3)

loglog(time,e2,'b','LineWidth',2)
loglog(time,errp_IAAFT_plus,'--b','LineWidth',3)
loglog(time,errp_IAAFT_minus,'--b','LineWidth',3)
lgd = legend('ORG', 'ORG ERRORBAR', 'ORG ERRORBAR','IAAFT','IAAFT ERRORBAR','IAAFT ERRORBAR', 'Location', 'northeast');
lgd.FontSize = 20;
%set(gca, 'XScale', 'log');
hold off
ax = gca;
ax.FontSize = 40;
grid on;
%ylim([0.18 0.85]);
%title('Heart Rate Multiscale Fuzzy Entropy');
xlabel('Time Scale [sec]');
ylabel('Fuzzy Entropy');
