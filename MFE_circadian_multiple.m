function [e_all, e_IAAFT_all] = MFE_circadian_multiple(data_cell, c, maxiter, m, factor, mf, rn, local, tau)
% MFE_circadian_multiple: Circadian Rhythm用の複数データを処理する関数

% 0405に実行したプログラム
% MFE_circadian_multiple(data_cell,10,50,2,10000,'Exponential',[0.2 2],0,1);
tic; % 時間計測用
% データの数
num_data = numel(data_cell);

% 出力用の変数を初期化
e_all = cell(1, num_data);
e_IAAFT_all = cell(1, num_data);
q = 0;

for data_index = 1
    
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
toc;    % 時間計測用
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


%%%%%%%%%%%
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


function e = fuzzymsentropy(input, m, mf, rn, local, tau, factor)
    y = input;
    y = y - mean(y);
    y = y / std(y);
    e = zeros(factor, 1);
    %pool = parpool(10);
    
    for i = 1:factor
        s = coarsegraining(y, i);
        sampe = FuzEn_MFs(s, m, mf, rn, local, tau);
        
        e(i) = sampe;
    end
    e = e';
    
end



% % Membership functions
% function c = Triangular(dist, rn)
%     c = zeros(size(dist));
%     c(dist <= rn) = 1 - dist(dist <= rn) ./ rn;
% end
% 
% function c = Trapezoidal(dist, rn)
%     c = zeros(size(dist));
%     c(dist <= rn) = 1;
%     c(dist <= 2 * rn & dist > rn) = 2 - dist(dist <= 2 * rn & dist > rn) ./ rn;
% end
% 
% function c = Z_shaped(dist, rn)
%     c = zeros(size(dist));
%     r1 = dist <= rn;
%     r2 = dist > rn & dist <= 1.5 * rn;
%     r3 = dist > 1.5 * rn & dist <= 2 * rn;
%     c(r1) = 1;
%     c(r2) = 1 - 2 .* ((dist(r2) - rn) ./ rn) .^ 2;
%     c(r3) = 2 .* ((dist(r3) - 2 * rn) ./ rn) .^ 2;
% end
% 
% function c = Bell_shaped(dist, rn)
%     c = 1 ./ (1 + abs(dist ./ rn(1)) .^ (2 * rn(2)));
% end
% 
% function c = Gaussian(dist, rn)
%     c = exp(-(dist ./ (sqrt(2) * rn)) .^ 2);
% end
% 
% function c = Constant_Gaussian(dist, rn)
%     c = ones(size(dist));
%     c(dist > rn) = exp(-log(2) .* ((dist(dist > rn) - rn) ./ rn) .^ 2);
% end
% 
% function c = Exponential(dist, rn)
%     c = exp(-dist .^ rn(2) ./ rn(1));
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function entr = FuzEn_MFs(ts, m, mf, rn, local, tau)
% 
% if nargin == 5, tau = 1; end
% if nargin == 4, local = 0; tau=1; end
% if nargin == 3, rn=0.2*std(ts);local = 0; tau=1; end
% 
% % parse inputs
% narginchk(6, 6);
% N = length(ts);
% 
% % normalization
% %ts = zscore(ts(:));
% 
% % reconstruction
% indm = hankel(1:N-m*tau, N-m*tau:N-tau);    % indexing elements for dim-m
% indm = indm(:, 1:tau:end);
% ym   = ts(indm);
% 
% inda = hankel(1:N-m*tau, N-m*tau:N);        % for dim-m+1
% inda = inda(:, 1:tau:end);
% ya   = ts(inda);
% 
% if local
%     ym = ym - mean(ym, 2)*ones(1, m);
%     ya = ya - mean(ya, 2)*ones(1, m+1);
% end
% 
% % inter-vector distance calculation in blocks
% blockSize = 1000; % Adjust block size based on available memory
% cheb_m = zeros(size(ym, 1), 1);
% cheb_a = zeros(size(ya, 1), 1);
% 
% for i = 1:blockSize:size(ym, 1)
%     endIndex = min(i+blockSize-1, size(ym, 1));
%     blockDistances = pdist2(ym(i:endIndex, :), ym, 'chebychev');
%     cheb_m(i:endIndex) = max(blockDistances, [], 2);
% end
% 
% for i = 1:blockSize:size(ya, 1)
%     endIndex = min(i+blockSize-1, size(ya, 1));
%     blockDistances = pdist2(ya(i:endIndex, :), ya, 'chebychev');
%     cheb_a(i:endIndex) = max(blockDistances, [], 2);
% end
% 
% cm = feval(mf, cheb_m, rn);
% ca = feval(mf, cheb_a, rn);
% 
% % output
% entr = -log(sum(ca) / sum(cm));
% clear indm ym inda ya cheb_m cheb_a cm ca;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 正しいやつ
% function entr = FuzEn_MFs(ts, m, mf, rn, local, tau)
% 
% if nargin == 5, tau = 1; end
% if nargin == 4, local = 0; tau=1; end
% if nargin == 3, rn=0.2*std(ts);local = 0; tau=1; end
% 
% % parse inputs
% narginchk(6, 6);
% N     = length(ts);
% 
% % normalization
% %ts = zscore(ts(:));
% 
% % reconstruction
% indm = hankel(1:N-m*tau, N-m*tau:N-tau);    % indexing elements for dim-m
% indm = indm(:, 1:tau:end);
% ym   = ts(indm);
% 
% inda = hankel(1:N-m*tau, N-m*tau:N);        % for dim-m+1
% inda = inda(:, 1:tau:end);
% ya   = ts(inda);
% 
% if local
%     ym = ym - mean(ym, 2)*ones(1, m);
%     ya = ya - mean(ya, 2)*ones(1, m+1);
% end
% 
% % inter-vector distance
% % if N < 1e4
% 
%     cheb = pdist(ym, 'chebychev'); % inf-norm
%     cm   = feval(mf, cheb, rn);
% 
% 
%     cheb = pdist(ya, 'chebychev');
%     ca   = feval(mf, cheb, rn);
% 
% % output
% entr = -log(sum(ca) / sum(cm));
% clear indm ym inda ya cheb cm ca;
% end
%%%%%%%%%%%%%
function entr = FuzEn_MFs(ts, m, mf, rn, local, tau)

if nargin == 5, tau = 1; end
if nargin == 4, local = 0; tau=1; end
if nargin == 3, rn=0.2*std(ts);local = 0; tau=1; end

% parse inputs
narginchk(6, 6);
N     = length(ts);

% normalization
%ts = zscore(ts(:));

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

% Inter-vector distance with batch processing
batch_size = 1000;  % Adjust this value based on available memory
num_batches = ceil(size(ym, 1) / batch_size);
cm = zeros(size(ym, 1), 1);
ca = zeros(size(ya, 1), 1);

for i = 1:num_batches
    batch_start = (i-1)*batch_size + 1;
    batch_end = min(i*batch_size, size(ym, 1));
    
    cheb_m_batch = pdist2(ym(batch_start:batch_end, :), ym, 'chebychev');
    cm(batch_start:batch_end) = max(cheb_m_batch, [], 2);
    
    cheb_a_batch = pdist2(ya(batch_start:batch_end, :), ya, 'chebychev');
    ca(batch_start:batch_end) = max(cheb_a_batch, [], 2);
end

cm = feval(mf, cm, rn);
ca = feval(mf, ca, rn);

% output
entr = -log(sum(ca) / sum(cm));
clear indm ym inda ya cheb cm ca;
end


%%%%%%%%%%%%
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

function c = Exponential(dist, rn)
    c = exp(-dist .^ rn(2) ./ rn(1));
end


