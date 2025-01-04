function [e_all, e_IAAFT_all] = MFE_circadian_multiple_1029_factor_cut(data_cell, c, maxiter, m, factor, mf, rn, local, tau, num)
    tic; % Start timer
    num_data = numel(data_cell); % Number of data cells
    e_all = cell(1, num_data);
    e_IAAFT_all = cell(1, num_data);

    for data_index = 1:1
        data = data_cell{data_index};
        data_l = length(data);
        disp('Start fuzzy')
        [e, e_IAAFT] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau, data_l, num);
        e_all{data_index} = e;
        e_IAAFT_all{data_index} = e_IAAFT;
        disp("Progress (15 levels):")
        disp(data_index)
    end

    toc; % End timer
end

function [e1, e_IAAFT] = MFE_circadian_single(data, c, maxiter, m, factor, mf, rn, local, tau, data_l, num)
    e1 = fuzzymsentropy(data, m, mf, rn, local, tau, factor, num);
    disp('Original completed');
    e2 = zeros(factor - num, c, 'gpuArray');
    [s, ~] = IAAFT(data, c, maxiter);
    disp('Surrogate completed');
    for i = 1:c
        e2(:, i) = fuzzymsentropy(s(:, i), m, mf, rn, local, tau, factor, num);
        disp("Fuzzy entropy progress for surrogate data")
        disp(i)
    end
    e2 = gather(e2)';
    e_IAAFT = e2;
    plot_MFE_graph(e1, e_IAAFT, data_l, num);
end

function plot_MFE_graph(e1, e2, data_l, num)
    time_length = data_l * 5;
    factor = size(e2, 2);
    time_s = zeros(factor, 1);
    time = zeros(factor, 1);
    for i = (num + 1):factor
        time_s(i) = data_l / i;
        time(i - num) = (time_length / time_s(i));
    end
    figure;
    plot(time, e1, 'r')
    hold on
    e2_mean = mean(e2, 1);
    e2_std = std(e2, 0, 1);
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

function e = fuzzymsentropy(input, m, mf, rn, local, tau, factor, num)
    y = input;
    y = y - mean(y);
    y = y / std(y);
    e = zeros(factor - num, 1, 'gpuArray');
    for i = (num + 1):factor
        s = coarsegraining(y, i);
        sampe = FuzEn_MFs(s, m, mf, rn, local, tau);
        e(i - num) = sampe;
    end
    e = gather(e)';
    disp('Progress')
end

function s = coarsegraining(inputSignal, scaleFactor)
    signalLength = length(inputSignal);
    newLength = floor(signalLength / scaleFactor);
    s = zeros(1, newLength, 'gpuArray');
    for i = 1:newLength
        startIndex = (i - 1) * scaleFactor + 1;
        endIndex = i * scaleFactor;
        s(i) = mean(inputSignal(startIndex:endIndex));
    end
end

function [s, r] = IAAFT(data, c, maxiter)
    data_length = length(data);
    s = gpuArray.zeros(data_length, c);
    for j = 1:c
        r = gpuArray(data(randperm(data_length)));
        A = abs(fft(data));
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

% function entr = FuzEn_MFs(ts, m, mf, rn, local, tau)
%     if nargin == 5, tau = 1; end
%     if nargin == 4, local = 0; tau=1; end
%     if nargin == 3, rn=0.2*std(ts);local = 0; tau=1; end
%     narginchk(6, 6);
%     N = length(ts);
%     indm = hankel(1:N-m*tau, N-m*tau:N-tau);
%     indm = indm(:, 1:tau:end);
%     ym = ts(indm);
%     inda = hankel(1:N-m*tau, N-m*tau:N);
%     inda = inda(:, 1:tau:end);
%     ya = ts(inda);
%     if local
%         ym = ym - mean(ym, 2)*ones(1, m);
%         ya = ya - mean(ya, 2)*ones(1, m+1);
%     end
%     cheb = pdist(ym, 'chebychev', 'Smallest', 1, 'UseParallel', true, 'CacheSize', 1024);
%     cm = feval(mf, cheb, rn);
%     cheb = pdist(ya, 'chebychev', 'Smallest', 1, 'UseParallel', true, 'CacheSize', 1024);
%     ca = feval(mf, cheb, rn);
%     entr = -log(sum(ca) / sum(cm));
% end
function entr = FuzEn_MFs(ts, m, mf, rn, local, tau)
    % Set default parameters
    if nargin == 5, tau = 1; end
    if nargin == 4, local = 0; tau = 1; end
    if nargin == 3, rn = 0.2 * std(ts); local = 0; tau = 1; end
    narginchk(6, 6);
    
    % Move input to GPU
    ts = gpuArray(ts); 
    N = length(ts);
    
    % Create indices for coarse graining
    indm = hankel(1:N-m*tau, N-m*tau:N-tau);
    indm = indm(:, 1:tau:end);
    ym = ts(indm);
    
    inda = hankel(1:N-m*tau, N-m*tau:N);
    inda = inda(:, 1:tau:end);
    ya = ts(inda);
    
    % Mean subtraction if local
    if local
        ym = ym - mean(ym, 2) * ones(1, m, 'gpuArray');
        ya = ya - mean(ya, 2) * ones(1, m + 1, 'gpuArray');
    end
    
    % Custom GPU-based Chebyshev distance calculation
    cheb_m = custom_chebyshev_distance(ym);
    cm = feval(mf, cheb_m, rn);
    
    cheb_a = custom_chebyshev_distance(ya);
    ca = feval(mf, cheb_a, rn);
    
    % Calculate entropy
    entr = -log(sum(ca) / sum(cm));
    
    % Gather the result from GPU
    entr = gather(entr);
end

function dist = custom_chebyshev_distance(data)
    % Custom implementation of Chebyshev distance on GPU
    n = size(data, 1);
    dist = zeros(n, n, 'gpuArray');
    for i = 1:n
        for j = i+1:n
            dist(i, j) = max(abs(data(i, :) - data(j, :)));
            dist(j, i) = dist(i, j); % Distance matrix is symmetric
        end
    end
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

% 大体このメンバーシップ関数を使用している
function c = Exponential(dist, rn)
    c = exp(-dist .^ rn(2) ./ rn(1));
end

