function [e_all, e_IAAFT_all] = MFE_circadian_multiple(data_cell, c, maxiter, m, factor, mf, rn, local, tau, num)
    % マルチスケールファジーエントロピーの計算 (複数データ)
    tic;
    num_data = numel(data_cell);
    e_all = cell(1, num_data);
    e_IAAFT_all = cell(1, num_data);

    for data_index = 1:num_data
        fprintf('Processing data index %d/%d...\n', data_index, num_data);
        data = data_cell{data_index};
        data_l = length(data);

        [e, e_IAAFT] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau, num);
        e_all{data_index} = e;
        e_IAAFT_all{data_index} = e_IAAFT;
    end

    toc;
end

function [e1, e_IAAFT] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau, num)
    % シングルデータ用 MFE 計算
    e1 = fuzzymsentropy(data, m, mf, rn, local, tau, factor, num);

    e2 = zeros(factor - num, c, 'gpuArray');
    [s, ~] = IAAFT(data, c, maxiter);

    for i = 1:c
        e2(:, i) = fuzzymsentropy(s(:, i), m, mf, rn, local, tau, factor, num);
    end

    e_IAAFT = gather(e2');
end

function e = fuzzymsentropy(data, m, mf, rn, local, tau, factor, num)
    % マルチスケールファジーエントロピー
    data = (data - mean(data)) / std(data);
    e = zeros(factor - num, 1);
    for i = (num + 1):factor
        s = coarsegraining(data, i);
        e(i - num) = FuzEn_MFs(s, m, mf, rn, local, tau);
    end
end

function s = coarsegraining(data, scaleFactor)
    % 粗視化処理
    signalLength = length(data);
    newLength = floor(signalLength / scaleFactor);
    s = arrayfun(@(i) mean(data((i - 1) * scaleFactor + 1:i * scaleFactor)), 1:newLength);
end

function [s, r] = IAAFT(data, c, maxiter)
    % IAAFT サロゲートデータ生成
    data_length = length(data);
    s = gpuArray.zeros(data_length, c);
    A = abs(fft(data));

    for j = 1:c
        r = data(randperm(data_length));
        for iter = 1:maxiter
            S = fft(r);
            S = A .* exp(1i * angle(S));
            r = real(ifft(S));
            [~, I] = sort(r);
            [~, J] = sort(data);
            r(I) = data(J);
        end
        s(:, j) = r;
    end
end

function entr = FuzEn_MFs(ts, m, mf, rn, local, tau)
    % ファジーエントロピー計算
    N = length(ts);
    
    % ヘンケル行列の作成
    ym = hankel(ts(1:N-m*tau), ts(N-m*tau:N-tau));  % ymは (N-m*tau)x(m) の行列
    ya = hankel(ts(1:N-m*tau), ts(N-m*tau:N));      % yaは (N-m*tau)x(m+1) の行列
    
    if local
        % ローカル平均を引く
        ym = ym - mean(ym, 2);
        ya = ya - mean(ya, 2);
    end
    
    % chebyshev 距離の計算を並列化
    cheb_m = zeros(size(ym, 1), size(ym, 2));  % chebyshev 距離の配列
    cheb_a = zeros(size(ya, 1), size(ya, 2));  % chebyshev 距離の配列
    
    % parfor で並列処理
    parfor i = 1:size(ym, 2)
        % ymとyaの次元を合わせるためにreshapeを調整
        cheb_m(:, i) = max(abs(ym(:, i) - ym), [], 2);  % m次元の距離計算
        cheb_a(:, i) = max(abs(ya(:, i) - ya), [], 2);  % (m+1)次元の距離計算
    end
    
    % メンバーシップ関数の適用
    cm = feval(mf, cheb_m(:), rn);
    ca = feval(mf, cheb_a(:), rn);
    
    % ファジーエントロピーの計算
    entr = -log(sum(ca) / sum(cm));
end




% 大体このメンバーシップ関数を使用している
function c = Exponential(dist, rn)
    c = exp(-dist .^ rn(2) ./ rn(1));
end

