function  pwelch_all_person_analysis(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10)


%rng default
fs=0.2;
t=0:1/fs:5-2/fs;

% [pxx,f] = pwelch(Z_all_ave,368377,184150,368377,fs);
%  sd_all=std(Z_all);
%  M_all=mean(Z_all);
% errp=M_all+sd_all;
% errm=M_all-sd_all;
% figure('Name','Allperson power','NumberTitle','off');
% %semilogy(f,pxx)
% loglog(f,M_all,'-k');
% hold on;
% loglog(f,errp,'--b');
% hold on;
% loglog(f,errm,'--b');
% grid on;
% ylim([100 1000000]);
% title('Allperson AllDay');
% xlabel('Frequency(Hz)');
% ylabel('Power');

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


[pxx1,f] = pwelch(Z_maru,122792,49116,122792,fs);
[pxx2,f] = pwelch(Z_hiki,122792,49116,122792,fs);
[pxx3,f] = pwelch(Z_yasu,122792,49116,122792,fs);
[pxx4,f] = pwelch(Z_kawa,122792,49116,122792,fs);
[pxx5,f] = pwelch(Z_shima,122792,49116,122792,fs);
[pxx6,f] = pwelch(Z_tuno,122792,49116,122792,fs);
[pxx7,f] = pwelch(Z_huka,122792,49116,122792,fs);
[pxx8,f] = pwelch(Z_huzi,122792,49116,122792,fs);
[pxx9,f] = pwelch(Z_mine,122792,49116,122792,fs);
[pxx10,f] = pwelch(Z_yama,122792,49116,122792,fs);

Z_all = (pxx1 + pxx2 + pxx3 + pxx4 + pxx5 + pxx6 + pxx7 + pxx8 + pxx9 + pxx10) / 10;

%[pxx,f] = pwelch(Z_all_ave,368377,184150,368377,fs);

 sd_all=std(Z_all);
 M_all=mean(Z_all);
errp=M_all+sd_all;
errm=M_all-sd_all;
figure('Name','Allperson power','NumberTitle','off');
%semilogy(f,pxx)
loglog(f,M_all,'-k');
hold on;
loglog(f,errp,'--b');
hold on;
loglog(f,errm,'--b');
grid on;
ylim([100 1000000]);
title('Allperson AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');
