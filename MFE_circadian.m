function [e1,e_IAAFT] = MFE_circadian(c,maxiter,m,factor,mf,rn,local,tau,data1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multiscale Fuzzy Entropy  専用の計算関数 (Circadian Rhythm用)


%IAAFT
% c: 元のデータに対して，位相をシャッフルしたデータの個数
% maxiter：データに対して，位相をシャッフルする回数


%fuzzymse(マルチスケールファジーエントロピー)
% m=2; %次元
% factor=3000;
% mf='Exponential';%指数関数(生体系はこれを使っていることが多い)
% rn=[0.2 2];%
% local=0;%global similarity
% tau=1;%時間遅延は一つ次にする(sampleEnと同じ),
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 使用する心拍変動データの初期化
y1 = data1;
l1 = length(y1);

%y1 = single(y1); % ここを変える可能性は高い(データを軽くしている)

% Fuzzy Entropy と IAAFT を行う上で必要な関数が格納されているディレクトリの指定 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath 'C:\Users\iw200\OneDrive\ドキュメント\MATLAB\Examples\R2023b\matlab\git_pupil'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IAAFTと元データのマルチスケールファジーエントロピー
input1 = y1;
e1 = fuzzymsentropy(input1,m,mf,rn,local,tau,factor);
e2 = zeros(factor,c);
[s,iter]=IAAFT(input1,c,maxiter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% サロゲートデータのマルチスケールファジーエントロピー
input2 = s;
for i = 1:c
%input2 = input2(:,i);
%e2(:,i) = msentropy_kai(input2(:,i),m,r,factor);
e2(:,i) = fuzzymsentropy(input2(:,i),m,mf,rn,local,tau,factor);

end
%転置行列(1000,10)→(10,1000) (factor,c)→(c,factor)
e2 = e2';
e_IAAFT = mean(e2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_length = l1 * 5;
for i=1:factor
time_s(i) = l1 / i;
time(i) = time_length / time_s(i);
end
% グラフの表示(マルチスケールファジーエントロピー)
%figure('Name','HeartRate of Multiscale Fuzzy Entropy','NumberTitle','off')
figure;
plot(time,e1,'r')
hold on
errorbar(time,mean(e2),std(e2),'b');
legend('ORG','IAAFT','Location','southeast');
set(gca, 'XScale', 'log');
hold off
title('Heart Rate Multiscale Fuzzy Entropy');

xlabel('Time Scale');
ylabel('Fuzzy Entropy');

