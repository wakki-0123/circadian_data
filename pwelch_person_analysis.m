function pwelch_person_analysis(data1,data2,data3,data4,data5,data6,data7,data8,data9,data10)

% 心拍変動のパワー解析用の関数(各被験者のパワー)

fs = 0.2; % sampling frequency
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

% 各被験者のパワー解析の結果の描画
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yasuda
figure('Name','yasuda power','NumberTitle','off');
loglog(f1,pxx1,'-k');
grid on;
ylim([100 1000000]);
title('yasuda AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% hikiti
figure('Name','hikiti power','NumberTitle','off');
loglog(f2,pxx2,'-k');
grid on;
ylim([100 1000000]);
title('hikiti AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tunoda
figure('Name','tunoda power','NumberTitle','off');
loglog(f3,pxx3,'-k');
grid on;
ylim([100 1000000]);
title('tunoda AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maruyama
figure('Name','maruyama power','NumberTitle','off');
loglog(f4,pxx4,'-k');
grid on;
ylim([100 1000000]);
title('maruyama AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% yamamoto
figure('Name','yamamoto power','NumberTitle','off');
loglog(f5,pxx5,'-k');
grid on;
ylim([100 1000000]);
title('yamamoto AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fukazawa
figure('Name','fukazawa power','NumberTitle','off');
loglog(f6,pxx6,'-k');
grid on;
ylim([100 1000000]);
title('fukazawa AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kawate
figure('Name','kawate power','NumberTitle','off');
loglog(f7,pxx7,'-k');
grid on;
ylim([100 1000000]);
title('kawate AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simada
figure('Name','simada power','NumberTitle','off');
loglog(f8,pxx8,'-k');
grid on;
ylim([100 1000000]);
title('simada AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fujimoto
figure('Name','fujimoto power','NumberTitle','off');
loglog(f9,pxx9,'-k');
grid on;
ylim([100 1000000]);
title('fujimoto AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minegisi
figure('Name','minegisi power','NumberTitle','off');
loglog(f10,pxx10,'-k');
grid on;
ylim([100 1000000]);
title('minegisi AllDay');
xlabel('Frequency(Hz)');
ylabel('Power');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


