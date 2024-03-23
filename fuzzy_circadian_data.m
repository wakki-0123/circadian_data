function y1 = fuzzy_circadian_data(c,maxiter,m,factor,mf,rn,local,tau)

% 心拍変動のplotとmultiscale fuzzy entropy のプログラム
% 元データの前処理とIAAFT、そしてファジーエントロピーを用いたマルチスケールファジーエントロピー解析

%前処理


%IAAFT
% c: 元のデータに対して，位相をシャッフルしたデータの個数
% maxiter：データに対して，位相をシャッフルする回数


%fuzzymse(マルチスケールファジーエントロピー)
% m=2; %次元
% factor=3000;
% mf='Exponential';%指数関数(生体系はこれを使っていることが多い)
% rn=[0.2 2];%
% local=0;%global similarity
% tau=1;%時間遅延は一つ次にする(sampleEnと同じ),'

file = 'fitbit_heartrate_修論用_安田.xlsx';
%データの読み込み

% Excelファイルから数値データを読み込む
[numData, ~, ~] = xlsread('fitbit_heartrate_修論用_安田.xlsx');

% データの確認
%disp(numData(:,3));





%前処理
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = numData(:,3);
%y = medfilt1(data,thredshold,'omitnan','truncate');%メジアン補間(NaNの処理)
%y = fillmissing(y,'linear');%線形補間
%y = filloutliers(y,'linear','movmedian',thredshold);



figure;
plot(y)

legend('ORG','Location','southeast');
title('original heartrate data');
xlabel('sample');
ylabel('heartrate');

l = length(y);
t_original = 1:1:l;
% 新しいサンプリング周波数を指定
new_sampling_frequency = 0.2;  % 新しいサンプリング周波数 (Hz)

% 新しい時間ベクトルを自動的に計算
new_time_vector = 1:1/new_sampling_frequency:max(t_original);

% resample関数を使用してデータをリサンプリング
y_resampled = interp1(t_original, y, new_time_vector);



y1 = y_resampled;
y1 = zscore(y1);


%y = data; %元のデータを表示したいときはコメントを外す(欠損値だらけ)
% 前処理後の時系列データのplot
figure;
plot(y1)

legend('ORG','Location','southeast');
title('resample heartrate data');
xlabel('time[s]');
ylabel('heartrate');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath 'C:\Users\iw200\Documents\MATLAB\eeg\pupil'
% IAAFTと元データのマルチスケールファジーエントロピー
input1 = y1;
%e1 = msentropy_kai(input1,m,r,factor);

%e = fuzzymsentropy(input,m,mf,rn,local,tau,factor)

e1 = fuzzymsentropy(input1,m,mf,rn,local,tau,factor);
e2 = zeros(factor,c);
[s,iter]=IAAFT(input1,c,maxiter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% input2 = y1;
% %e1 = msentropy_kai(input1,m,r,factor);
% 
% %e = fuzzymsentropy(input,m,mf,rn,local,tau,factor)
% 
% e3 = fuzzymsentropy(input2,m,mf,rn,local,tau,factor);
% e4 = zeros(factor,c);
% [s1,iter1]=IAAFT(input2,c,maxiter);

% figure;
% plot(s);
% title('IAAFT')
% xlabel('sample')
% ylabel('pipil')
% legend('surrogate signal','Location','southeast')
% disp(size(s));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% サロゲートデータのマルチスケールファジーエントロピー
input2 = s;
for i = 1:c
%input2 = input2(:,i);
%e2(:,i) = msentropy_kai(input2(:,i),m,r,factor);
e2(:,i) = fuzzymsentropy(input2(:,i),m,mf,rn,local,tau,factor);

end
%転置行列(1000,10)→(10,1000) (factor,c)→(c,factor)
e2 = e2';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% サロゲートデータのマルチスケールファジーエントロピー
% input3 = s1;
% for i = 1:c
% %input2 = input2(:,i);
% %e2(:,i) = msentropy_kai(input2(:,i),m,r,factor);
% e4(:,i) = fuzzymsentropy(input3(:,i),m,mf,rn,local,tau,factor);

% end
% %転置行列(1000,10)→(10,1000) (factor,c)→(c,factor)
% e4 = e4';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% グラフの表示(マルチスケールファジーエントロピー)
figure(3)
plot(1:1:factor,e1,'r')
hold on
errorbar(1:1:factor,mean(e2),std(e2),'b');
legend('ORG','IAAFT','Location','southeast');
set(gca, 'XScale', 'log');
hold off
title('left');

xlabel('time scale[sec]');
ylabel('fuzzy entropy');