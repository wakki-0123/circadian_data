function [I_Cx, I_Cxx, I_xx, H_x, H_xx, H_C] = Fuzzy_MI(data)
% Optimized Fuzzy Mutual Information and Fuzzy Entropy calculation

data(:, end) = grp2idx(data(:, end));
NC = max(data(:, end));
NF = size(data, 2) - 1;
NP = size(data, 1);

R = NC;
data1 = mapstd(data(:, 1:NF)', 0, 1);
data(:, 1:NF) = data1';
clear data1

H_x = zeros(1, NF);
H_xx = zeros(NF, NF);
I_xx = zeros(NF, NF);
I_Cx = zeros(1, NF);
I_Cxx = zeros(NF, NF);
U = zeros(NC, NP, NF);
Pr_C = histcounts(data(:, end), NC) / NP;

H_C = -sum(Pr_C .* log(Pr_C + eps));    % entropy of class label

C = MeanClust(data);            % Mean of the different classes
expo = 1.152;                    % fuzzifier
for i = 1:NF
    distA = ((Rdistfcm(C(:, i), data(:, i)) + eps));
    U_new = ((distA) .^ (-2 / (expo - 1)));
    U_new = U_new ./ sum(U_new, 1);
    U(:, :, i) = U_new;
    temp = (U_new * U_new') / NP;
    temp = temp(temp > 0);
    H_x(i) = -sum(temp .* log(temp));
    Ps_Cx = zeros(NC, NC, NC);
    for j = 1:NC
        matC = data(:, end) == j;
        Ps_Cx(:, :, j) = ((U(:, :, i)' .* matC')' * (U(:, :, i)' .* matC')) / NP;
    end
    temp = Ps_Cx(Ps_Cx > 0);
    H_Cx = -sum(temp .* log(temp));
    I_Cx(i) = H_x(i) + H_C - H_Cx;
end

tempC = zeros(NC, R, R);

for n = 1:NF-1
    for m = n+1:NF
        xxa = double(squeeze(U(:, :, n)));
        xxb = double(squeeze(U(:, :, m)));
        temp = (xxa * xxb') / NP;
        temp = temp(temp > 0);
        H_xx(n, m) = -sum(temp .* log(temp));
        I_xx(n, m) = H_x(n) + H_x(m) - H_xx(n, m);
        for k = 1:NC
            matC = data(:, end) == k;
            tempC(k, :, :) = ((xxa' .* matC')' * (xxb' .* matC')) / NP;
        end
        temp = tempC(tempC > 0);
        H_Cxx = -sum(temp .* log(temp));
        I_Cxx(n, m) = H_xx(n, m) + H_C - H_Cxx;
    end
end

end

function center = MeanClust(data)
M = max(data(:, end));
center = zeros(M, size(data, 2) - 1);
for i = 1:M
    classs = data(data(:, end) == i, 1:end-1);
    if size(classs, 1) > 1
        center(i, :) = mean(classs);
    else
        center(i, :) = classs;
    end
end
end

function out = Rdistfcm(center, data)
out = zeros(size(center, 1), size(data, 1));
if size(center, 2) > 1
    for k = 1:size(center, 1)
        out(k, :) = sqrt(sum((((data - center(k, :)).^2)' ./ var(data))'));
    end
else
    for k = 1:size(center, 1)
        out(k, :) = abs((center(k) - data) ./ var(data))';
    end
end
end
