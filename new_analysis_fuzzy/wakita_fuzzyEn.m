function [e_all, e_IAAFT_all] = wakita_fuzzyEn(data_cell, c, maxiter, m, factor, rn, local, tau)
    % Function to process multiple data for Circadian Rhythm using Multiscale Fuzzy Entropy

    tic; % Start timing

    % Number of data sets
    num_data = numel(data_cell);

    % Initialize output variables
    e_all = cell(1, num_data);
    e_IAAFT_all = cell(1, num_data);

    for data_index = 1:num_data
        % Get each data set
        data = data_cell{data_index};
        data_l = length(data);

        % Calculate Multiscale Fuzzy Entropy and IAAFT
        [e, e_IAAFT] = MFE_circadian_single(data, c, maxiter, m, factor, rn, local, tau, data_l);

        % Save results
        e_all{data_index} = e;
        e_IAAFT_all{data_index} = e_IAAFT;
    end

    toc; % End timing
end

function [e1, e_IAAFT] = MFE_circadian_single(data, c, maxiter, m, factor, rn, local, tau, data_l)
    % Function to calculate Multiscale Fuzzy Entropy and IAAFT for a single data set

    % Set the number of scales
    num = 7;

    % Calculate Multiscale Fuzzy Entropy
    e1 = fuzzymsentropy(data, m, rn, local, tau, factor, num);

    % Initialize for IAAFT surrogate data
    e2 = zeros(factor - num, c);
    [s, ~] = IAAFT(data, c, maxiter);

    % Calculate Multiscale Fuzzy Entropy for surrogate data
    for i = 1:c
        e_temp = fuzzymsentropy(s(:, i), m, rn, local, tau, factor, num);
        e2(:, i) = e_temp';
    end

    e_IAAFT = mean(e2, 2); % Take the mean across the surrogates

    % Plot the results
    plot_MFE_graph(e1, e2, data_l, num);
end

function plot_MFE_graph(e1, e2, data_l, num)
    % Function to plot Multiscale Fuzzy Entropy graph

    time_length = data_l * 5; % Total duration in seconds
    factor = size(e2, 1) + num;

    % Calculate time scales (excluding the first `num` scales)
    time_s = zeros(1, factor);
    time = zeros(1, factor - num);
    for i = (num + 1):factor
        time_s(i) = data_l / i; % Total number of samples
        time(i - num) = time_length / time_s(i); % Time scale
    end

    % Plot the graph
    figure;
    plot(time, e1, 'r')
    hold on
    errorbar(time, mean(e2, 2), std(e2, 0, 2), 'b');

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

function e = fuzzymsentropy(input, m, rn, local, tau, factor, num)
    % Function to calculate Multiscale Fuzzy Entropy

    addpath('./nobukawa_fuzzy_entropy/EntropyHub_v2.0.0/');
    y = input;
    y = y - mean(y);
    y = y / std(y);
    e = zeros(factor - num, 1);

    for i = (num + 1):factor
        s = coarsegraining(y, i);
        Mobj = MSobject('FuzzEn');
        e_temp = cMSEn(s, Mobj, 'Scales', factor);
        e(i - num, 1) = mean(e_temp); % Ensure compatibility
    end
    e = e';
end

function s = coarsegraining(inputSignal, scaleFactor)
    % Function to perform coarse-graining

    % Get the length of the input signal
    signalLength = length(inputSignal);

    % Calculate the number of new samples based on the scale factor
    newLength = floor(signalLength / scaleFactor);

    % Initialize the coarse-grained signal array
    s = zeros(1, newLength);

    % Perform coarse-graining
    for i = 1:newLength
        % Calculate each coarse-grained data point
        startIndex = (i - 1) * scaleFactor + 1;
        endIndex = i * scaleFactor;
        s(i) = mean(inputSignal(startIndex:endIndex));
    end
end

function [s, r] = IAAFT(data, c, maxiter)
    % Function to perform IAAFT surrogate analysis

    % Get the length of the data
    data_length = length(data);

    % Initialize the surrogate data array on GPU
    s = gpuArray.zeros(data_length, c);

    % Generate `c` surrogate data sets
    for j = 1:c
        % Shuffle the data randomly
        r = data(randperm(data_length));

        % Calculate the amplitude spectrum of the original data
        A = abs(fft(data));  % Use fft with gpuArray as input

        for iter = 1:maxiter
            % Perform Fourier transform on the surrogate data
            S = fft(r);

            % Replace amplitude with that of the original data
            S = A .* exp(1i * angle(S));

            % Perform inverse Fourier transform to return to time domain
            r = real(ifft(S));

            % Rank ordering
            [~, I] = sort(r);
            [~, J] = sort(data);
            r(I) = data(J);
        end

        % Save the surrogate data
        s(:, j) = r;
    end
end

