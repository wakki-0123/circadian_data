function pwelch_wakita_analysis(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10)

% 心拍変動のパワー解析用の関数(全被験者の平均のパワー)

fs=0.2; % sampling frequency
% 0323現在のワークスペースの変数に従って書いています
Z_yasu = data1;
Z_hiki = data2;
Z_tuno = data3;
Z_maru = data4;
Z_yama = data5;
Z_huka = data6;
Z_kawa = data7;
Z_shima = data8;
Z_huzi = data9;
Z_mine = data10;

[pxx1,f1] = pwelch(Z_yasu,[],[],[],fs);
[pxx2,f2] = pwelch(Z_hiki,[],[],[],fs);
[pxx3,f3] = pwelch(Z_tuno,[],[],[],fs);
[pxx4,f4] = pwelch(Z_maru,[],[],[],fs);
[pxx5,f5] = pwelch(Z_yama,[],[],[],fs);
[pxx6,f6] = pwelch(Z_huka,[],[],[],fs);
[pxx7,f7] = pwelch(Z_kawa,[],[],[],fs);
[pxx8,f8] = pwelch(Z_shima,[],[],[],fs);
[pxx9,f9] = pwelch(Z_huzi,[],[],[],fs);
[pxx10,f10] = pwelch(Z_mine,[],[],[],fs);

pxx_all = (pxx1 + pxx2 + pxx3 + pxx4 + pxx5 + pxx6 + pxx7 + pxx8 + pxx9 + pxx10) / 10; % average of power
f_all = (f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10) / 10; % average of frequency


sd_all=std(pxx_all); % 標準偏差
%pxx_all=mean(pxx_all);
errp=pxx_all+sd_all; % エラーバー(上)
errm=pxx_all-sd_all; % エラーバー(下)
figure('Name','Allperson power','NumberTitle','off');
%semilogy(f,pxx)
loglog(f_all,pxx_all,'-k','LineWidth',4);
hold on;
loglog(f_all,errp,'--b','LineWidth',2);
hold on;
loglog(f_all,errm,'--b','LineWidth',2);
ax = gca;
ax.FontSize = 40;


grid on;
ylim([100 1000000]);

lgd = legend('MEAN','STD','Location','southeast');
lgd.FontSize = 40;
xlabel('Frequency(Hz)');
ylabel('Power');

ax = gca;
ax.FontSize = 40;
