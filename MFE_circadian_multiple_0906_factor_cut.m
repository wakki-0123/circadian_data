function [e_all, e_IAAFT_all] = MFE_circadian_multiple_0906_factor_cut(data_cell, c, maxiter, m, factor, mf, rn, local, tau,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% プログラムの説明
% MFE_circadian_multiple: Circadian Rhythm用の複数データを処理する関数
% 0405に実行したプログラム
% MFE_circadian_multiple(data_cell,10,50,2,10000,'Exponential',[0.2 2],0,1);
tic; % 時間計測用
% データの数
% num = 7 35sec cut　0918時点
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_data = numel(data_cell); % cellデータの数
% 出力用の変数を初期化 % 1 ← num_data(右の)
e_all = cell(1, num_data);
e_IAAFT_all = cell(1, num_data);
q = 0; % カウンタ変数

for data_index = 1:num_data
    % 各データを取得
    %q = q + 1;
    data = data_cell{data_index};
    data_l = length(data);
    % マルチスケールファジーエントロピーとIAAFTを計算
    disp('start fuzzy')
    [e, e_IAAFT] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau, data_l,num);
    % 結果を保存
    e_all{data_index} = e;
    e_IAAFT_all{data_index} = e_IAAFT;
    disp("現在の稼働状況15段階 (全体でどのくらい進んでいるか)")
    %disp(q)
end
%delete(pool);
toc;    % 時間計測用
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [e1, e_IAAFT] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau, data_l,num)
% MFE_circadian_single: 単一のデータに対してマルチスケールファジーエントロピーとIAAFTを計算する関数
% マルチスケールファジーエントロピーの計算
e1 = fuzzymsentropy(data, m, mf, rn, local, tau, factor,num);
disp('オリジナル完了');
% IAAFTを実行し、その結果のマルチスケールファジーエントロピーを計算
% 予備的な変数設定
%num = 7;
e2 = zeros(factor - num, c); % 修正: e2 のサイズを factor - num x c に設定
[s, ~] = IAAFT(data, c, maxiter);
disp('サロゲート完了');
% サロゲートデータのマルチスケールファジーエントロピーを計算
for i = 1:c
    e2(:, i) = fuzzymsentropy(s(:, i), m, mf, rn, local, tau, factor,num);
    disp("現在のサロゲートデータに関するファジーエントロピーの様子")
    disp(i)
end
% 転置行列
e2 = e2';
e_IAAFT = e2;
% グラフの表示
plot_MFE_graph(e1, e_IAAFT, data_l,num);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_MFE_graph(e1, e2, data_l,num)
    % マルチスケールファジーエントロピーのグラフをプロットする関数
    time_length = data_l * 5; % 全部の区間の秒数
    factor = size(e2, 2); % e2 の列数を取得
    %num = 7; % factor cut
    % 時間スケールを計算(num以前は除く)
    time_s = zeros(factor, 1);
    time = zeros(factor, 1);
    for i = (num + 1):factor
        time_s(i) = data_l / i; % 合計サンプルの個数
        time(i - num) = (time_length / time_s(i)); % タイムスケール
    end
    % グラフの表示
    figure;
    plot(time, e1, 'r')
    hold on
    % e2 の平均と標準偏差を計算
    e2_mean = mean(e2, 1); % 各列の平均
    e2_std = std(e2, 0, 1); % 各列の標準偏差
    % errorbar のプロット
    errorbar(time, e2_mean, e2_std, 'b');
    lgd = legend('ORG', 'IAAFT', 'Location', 'southeast');
    lgd.FontSize = 40;
    set(gca, 'XScale', 'log');
    ax = gca;
    ax.FontSize = 40;
    hold off
    title('Heart Rate Multiscale Fuzzy Entropy');
    xlabel('Time Scale');
    ylabel('Fuzzy Entropy');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function e = fuzzymsentropy(input, m, mf, rn, local, tau, factor,num)
    y = input;
    y = y - mean(y);
    y = y / std(y);
    % num は計算に使用しない部分を除外するための閾値
    %num = 7; % factor cut
    % e のサイズを factor - num に設定
    e = zeros(factor - num, 1);
    for i = (num + 1):factor
        s = coarsegraining(y, i);
        sampe = FuzEn_MFs(s, m, mf, rn, local, tau);
        % インデックスの調整
        e(i - num) = sampe;
    end
    e = e';
    disp('進捗状況')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 正しいやつ、fuzzy entropy計算の根幹部分
function entr = FuzEn_MFs(ts, m, mf, rn, local, tau)
if nargin == 5, tau = 1; end
if nargin == 4, local = 0; tau=1; end
if nargin == 3, rn=0.2*std(ts);local = 0; tau=1; end
% parse inputs
narginchk(6, 6);
N     = length(ts);
% reconstruction
indm = hankel(1:N-m*tau, N-m*tau:N-tau);    % indexing elements for dim-m
indm = indm(:, 1:tau:end);
ym   = ts(indm);
inda = hankel(1:N-m*tau, N-m*tau:N);        % for dim-m+1
inda = inda(:, 1:tau:end);
ya   = ts(inda);
if local
    ym = ym - mean(ym, 2)*ones(1, m);
    ya = ya - mean(ya, 2)*ones(1, m+1);
end
% inter-vector distance
% if N < 1e4
    cheb = pdist(ym, 'chebychev','CacheSize',1024); % inf-norm
    cm   = feval(mf, cheb, rn);
    cheb = pdist(ya, 'chebychev','CacheSize',1024);
    ca   = feval(mf, cheb, rn);
% output
entr = -log(sum(ca) / sum(cm));
clear indm ym inda ya cheb cm ca;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%membership functions
 function c = Triangular(dist, rn)
     c = zeros(size(dist));
     c(dist <= rn) = 1 - dist(dist <= rn) ./ rn;
 end
 function c = Trapezoidal(dist, rn)
    c = zeros(size(dist));
    c(dist <= rn) = 1;
    c(dist <= 2 * rn & dist > rn) = 2 - dist(dist <= 2 * rn & dist > rn) ./ rn;
end
function c = Z_shaped(dist, rn)
    c = zeros(size(dist));
    r1 = dist <= rn;
    r2 = dist > rn & dist <= 1.5 * rn;
    r3 = dist > 1.5 * rn & dist <= 2 * rn;
    c(r1) = 1;
    c(r2) = 1 - 2 .* ((dist(r2) - rn) ./ rn) .^ 2;
    c(r3) = 2 .* ((dist(r3) - 2 * rn) ./ rn) .^ 2;
end
function c = Bell_shaped(dist, rn)
    c = 1 ./ (1 + abs(dist ./ rn(1)) .^ (2 * rn(2)));
end
function c = Gaussian(dist, rn)
    c = exp(-(dist ./ (sqrt(2) * rn)) .^ 2);
end
function c = Constant_Gaussian(dist, rn)
    c = ones(size(dist));
    c(dist > rn) = exp(-log(2) .* ((dist(dist > rn) - rn) ./ rn) .^ 2);
end

% 大体このメンバーシップ関数を使用している
function c = Exponential(dist, rn)
    c = exp(-dist .^ rn(2) ./ rn(1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 粗視化処理
function s = coarsegraining(inputSignal, scaleFactor)
    % 入力シグナルの長さを取得
    signalLength = length(inputSignal);

    % スケールファクターに基づいて新しいサンプル数を計算
    newLength = floor(signalLength / scaleFactor);

    % 粗視化されたシグナルを格納する配列を初期化
    s = zeros(1, newLength);

    % コースグレーニング処理を実行
    for i = 1:newLength
        % 各粗視化されたデータポイントを計算
        startIndex = (i - 1) * scaleFactor + 1;
        endIndex = i * scaleFactor;
        s(i) = mean(inputSignal(startIndex:endIndex));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IAAFTサロゲート解析
function [s, r] = IAAFT(data, c, maxiter)
    % IAAFT: Iteratively Adjusted Amplitude and Frequency Transform
    
    % データの長さを取得
    data_length = length(data);
    
    % サロゲートデータ配列をGPU上で初期化
    s = gpuArray.zeros(data_length, c);
    
    % c個のサロゲートデータセットを生成
    for j = 1:c
        % データをランダムにシャッフル
        r = data(randperm(data_length));
        
        % オリジナルデータの振幅スペクトルを計算
        A = abs(fft(data));  % fft関数を使用してgpuArrayを引数にする
        
        for iter = 1:maxiter
            % サロゲートデータをフーリエ変換
            S = fft(r);
            
            % 振幅をオリジナルデータのものに置き換え
            S = A .* exp(1i * angle(S));
            
            % 逆フーリエ変換して時間領域に戻す
            r = real(ifft(S));
            
            % 順位付け
            [~, I] = sort(r);
            [~, J] = sort(data);
            r(I) = data(J);
        end
        
        % サロゲートデータを保存
        s(:, j) = r;
    end
end








