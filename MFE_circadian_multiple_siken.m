function [e_all, e_IAAFT_all] = MFE_circadian_multiple_siken(data_cell, c, maxiter, m, factor, mf, rn, local, tau)
    % MFE_circadian_multiple: Circadian Rhythm用の複数データを処理する関数

    tic

    % スライス変数としてカウンタを初期化
    
    % q_lock = parallel.pool.DataQueue;
    % % コールバック関数でカウンタをインクリメント
    % afterEach(q_lock, @(data) assignin('base', 'q_lock', evalin('base', 'q_lock') + data));

    % 結果を保存するためのプレースホルダ
    %numData = numel(data_cell);
    numData = 1;
    e_all = cell(1, numData);
    e_IAAFT_all = cell(1, numData);
    e_IAAFT_std = cell(1, numData);
    pool = parpool(2);


   
        parfor data_index = 1:2
        % 各データを取得
        data = data_cell{data_index};
        
        %data_l = length(data);
         

        % マルチスケールファジーエントロピーとIAAFTを計算
        [e, e_IAAFT_mean, e_std] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau);

        % 結果を保存
        e_all{data_index} = e; % original
        e_IAAFT_all{data_index} = e_IAAFT_mean; %IAAFT_mean
        e_IAAFT_std{data_index} = e_std; %IAAFT_std
        disp('進捗');
        disp(data_index);
        end
        % プールを閉じる
        delete(pool);
        
        

        % カウンタを更新
        % disp(q_lock)
        % send(q_lock, 1);
        
    
    
    
    data_l2 = factor; % 長さを求める
    disp(data_l2);


    ex = cell2mat(e_all);
    ex = reshape(ex, numData, data_l2);
    disp(size(ex));

    ey = cell2mat(e_IAAFT_std);
    ey = reshape(ey, numData, data_l2);
    disp(size(ey));


    ez = cell2mat(e_IAAFT_all);
    ez = reshape(ez, numData, data_l2);
    disp(size(ez));

    for j=1:2
    plot_MFE_graph(ex(j,:), ey(j,:), ez(j,:), data_l2);
    disp('表示回数')
    disp(j)
    end


    % カウンタの値を表示（並列処理後）
    %disp(q)

    toc

end


function [e1, e_IAAFT_mean, e_std] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau)
    % Transfer data to GPU
    data_gpu = gpuArray(data);

    % Calculate MFE using GPU
    e1 = fuzzymsentropy(data_gpu, m, mf, rn, local, tau, factor);

    % Perform IAAFT and calculate MFE on surrogates using GPU
    e2 = zeros(factor, c, 'gpuArray');
    [s, ~] = IAAFT(data_gpu, c, maxiter);

    for i = 1:c
        e2(:, i) = fuzzymsentropy(s(:, i), m, mf, rn, local, tau, factor);
        disp("現在のサロゲートデータに関するファジーエントロピーの様子")
        disp(i)
    end

    % Compute mean and standard deviation on the GPU
    e2 = gather(e2'); % gather to bring it back to CPU if needed
    e_IAAFT_mean = mean(e2, 1);
    e_std = std(e2, 0, 1);
end

function plot_MFE_graph(e1, e2, e3, data_l)
% e1:original e2:std, e3:mean
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
errorbar(time, e3, e2, 'b');
legend('ORG', 'IAAFT', 'Location', 'southeast');
set(gca, 'XScale', 'log');
hold off
title('Heart Rate Multiscale Fuzzy Entropy');
xlabel('Time Scale');
ylabel('Fuzzy Entropy');
%ylim([0 2])

end

function e = fuzzymsentropy(input, m, mf, rn, local, tau, factor)
    y = gpuArray(input); % Move input data to GPU
    y = y - mean(y);
    y = y / std(y);
    e = zeros(factor, 1, 'gpuArray'); % Use GPU for storage
    
    for i = 1:factor
        s = coarsegraining(y, i); % Modify coarsegraining if needed for GPU
        sampe = FuzEn_MFs(s, m, mf, rn, local, tau);
        e(i) = sampe;
    end
    e = gather(e'); % gather results back to CPU if needed
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%正しいやつ
function entr = FuzEn_MFs(ts, m, mf, rn, local, tau)
    if nargin == 5, tau = 1; end
    if nargin == 4, local = 0; tau=1; end
    if nargin == 3, rn=0.2*std(ts); local = 0; tau=1; end

    % Normalization
    ts_gpu = gpuArray(ts);  % Move to GPU

    % Reconstruction
    indm = hankel(1:length(ts_gpu)-m*tau, length(ts_gpu)-m*tau:length(ts_gpu)-tau); % indexing elements for dim-m
    indm = gpuArray(indm(:, 1:tau:end));
    ym = ts_gpu(indm);

    inda = hankel(1:length(ts_gpu)-m*tau, length(ts_gpu)-m*tau:length(ts_gpu)); % for dim-m+1
    inda = gpuArray(inda(:, 1:tau:end));
    ya = ts_gpu(inda);

    if local
        ym = ym - mean(ym, 2)*ones(1, m, 'gpuArray');
        ya = ya - mean(ya, 2)*ones(1, m+1, 'gpuArray');
    end

    % Inter-vector distance (on GPU)
    cheb_m = pdist(ym, 'chebychev', 'gpuArray'); % inf-norm
    cm = feval(mf, cheb_m, rn);

    cheb_a = pdist(ya, 'chebychev', 'gpuArray');
    ca = feval(mf, cheb_a, rn);

    % Output
    entr = -log(sum(ca) / sum(cm));
    entr = gather(entr); % gather results back to CPU if needed
end



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
% % Inter-vector distance with batch processing
% batch_size = 1000;  % Adjust this value based on available memory
% num_batches = ceil(size(ym, 1) / batch_size);
% cm = zeros(size(ym, 1), 1);
% ca = zeros(size(ya, 1), 1);
% 
% for i = 1:num_batches
%     batch_start = (i-1)*batch_size + 1;
%     batch_end = min(i*batch_size, size(ym, 1));
% 
%     cheb_m_batch = pdist2(ym(batch_start:batch_end, :), ym, 'chebychev');
%     cm(batch_start:batch_end) = max(cheb_m_batch, [], 2);
% 
%     cheb_a_batch = pdist2(ya(batch_start:batch_end, :), ya, 'chebychev');
%     ca(batch_start:batch_end) = max(cheb_a_batch, [], 2);
% end
% 
% cm = feval(mf, cm, rn);
% ca = feval(mf, ca, rn);
% 
% % output
% entr = -log(sum(ca) / sum(cm));
% clear indm ym inda ya cheb cm ca;
% end

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
