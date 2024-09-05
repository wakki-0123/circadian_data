function [e_all, e_IAAFT_all] = MFE_circadian_multiple_siken2(data_cell, c, maxiter, m, factor, mf, rn, local, tau)
    % MFE_circadian_multiple: Circadian Rhythm用の複数データを処理する関数

    tic
    numData = numel(data_cell);
    e_all = cell(1, numData);
    e_IAAFT_all = cell(1, numData);
    e_IAAFT_std = cell(1, numData);
    
    pool = parpool('local', 2);  % 並列処理用のプールを作成
    %progress_bar = waitbar(0, '計算中...');  % プログレスバーの表示
    
    parfor data_index = 1:numData
        % 各データを取得
        data = data_cell{data_index};
        
        % マルチスケールファジーエントロピーとIAAFTを計算
        [e, e_IAAFT_mean, e_std] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau);

        % 結果を保存
        e_all{data_index} = e; % original
        e_IAAFT_all{data_index} = e_IAAFT_mean; %IAAFT_mean
        e_IAAFT_std{data_index} = e_std; %IAAFT_std
        
        % 進捗表示
        fprintf('データセット %d/%d 処理完了\n', data_index, numData);
        
    end

    delete(pool);  % プールを閉じる
    %close(progress_bar);  % プログレスバーを閉じる
    
    % ここから後は計算結果をプロットする処理です
    data_l2 = factor;
    ex = cell2mat(e_all);
    ex = reshape(ex, numData, data_l2);
    
    ey = cell2mat(e_IAAFT_std);
    ey = reshape(ey, numData, data_l2);
    
    ez = cell2mat(e_IAAFT_all);
    ez = reshape(ez, numData, data_l2);

    for j = 1:numData
        plot_MFE_graph(ex(j,:), ey(j,:), ez(j,:), data_l2);
        fprintf('グラフ表示 %d/%d\n', j, numData);
    end

    toc
end

function [e1, e_IAAFT_mean, e_std] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau)
    data_gpu = gpuArray(data);  % GPU上でデータを処理

    e1 = fuzzymsentropy(data_gpu, m, mf, rn, local, tau, factor);

    e2 = zeros(factor, c, 'gpuArray');
    [s, ~] = IAAFT(data_gpu, c, maxiter);

    for i = 1:c
        e2(:, i) = fuzzymsentropy(s(:, i), m, mf, rn, local, tau, factor);
        fprintf('サロゲートデータ %d/%d 処理中...\n', i, c);  % サロゲートデータの進捗表示
    end

    e2 = gather(e2');
    e_IAAFT_mean = mean(e2, 1);
    e_std = std(e2, 0, 1);
end

function plot_MFE_graph(e1, e2, e3, data_l)
    time_length = data_l * 5;
    factor = size(e2, 2);
    time = time_length ./ (data_l ./ (1:factor));

    figure;
    plot(time, e1, 'r')
    hold on
    errorbar(time, e3, e2, 'b');
    legend('ORG', 'IAAFT', 'Location', 'southeast');
    set(gca, 'XScale', 'log');
    hold off
    title('Heart Rate Multiscale Fuzzy Entropy');
    xlabel('Time Scale');
    ylabel('Fuzzy Entropy');
end

function e = fuzzymsentropy(input, m, mf, rn, local, tau, factor)
    y = gpuArray(input);
    y = y - mean(y);
    y = y / std(y);
    e = zeros(factor, 1, 'gpuArray');
    
    for i = 1:factor
        s = coarsegraining(y, i);
        sampe = FuzEn_MFs(s, m, mf, rn, local, tau);
        e(i) = sampe;
    end
    e = gather(e');
end

% function entr = FuzEn_MFs(ts, m, mf, rn, local, tau)
%     if nargin == 5, tau = 1; end
%     if nargin == 4, local = 0; tau=1; end
%     if nargin == 3, rn = 0.2 * std(ts); local = 0; tau = 1; end
% 
%     ts_gpu = gpuArray(ts); % Move to GPU
%     N = length(ts_gpu);
% 
%     indm = hankel(1:N-m*tau, N-m*tau:N-tau);
%     indm = gpuArray(indm(:, 1:tau:end));
%     ym = ts_gpu(indm);
% 
%     inda = hankel(1:N-m*tau, N-m*tau:N);
%     inda = gpuArray(inda(:, 1:tau:end));
%     ya = ts_gpu(inda);
% 
%     if local
%         ym = ym - mean(ym, 2) * ones(1, m, 'gpuArray');
%         ya = ya - mean(ya, 2) * ones(1, m + 1, 'gpuArray');
%     end
% 
%     cheb_m = pdist2(ym, ym, 'chebychev', 'gpuArray');
%     cm = feval(mf, cheb_m, rn);
% 
%     cheb_a = pdist2(ya, ya, 'chebychev', 'gpuArray');
%     ca = feval(mf, cheb_a, rn);
% 
%     entr = -log(sum(ca) / sum(cm));
%     entr = gather(entr); % gather results back to CPU if needed
% end


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
