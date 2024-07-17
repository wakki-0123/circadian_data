function entr = FuzEn_MFs(ts, m, mf, rn, local, tau)
% This function calculates fuzzy entropy (FuzEn) of a univariate
% signal ts, using different fuzzy membership functions (MFs)

% Inputs:
%   ts      - time-series - a vector of size 1 x N (the number of sample points)
%   m       - embedding dimension
%   mf      - membership function, e.g., 'Triangular', 'Trapezoidal', etc.
%   rn      - threshold r and order n (scalar or vector based on mf)
%             scalar: threshold
%             vector: [threshold r, order n]
%   local   - local similarity (1) or global similarity (0)
%   tau     - time delay

% Reference:
%   [1] H. Azami, P. Li, S. Arnold, J. Escudero, and A. Humeau-Heurtier,
%       "Fuzzy Entropy Metrics for the Analysis of Biomedical Signals:
%       Assessment and Comparison", IEEE ACCESS, 2019.

% Authors: Peng Li and Hamed Azami
% Emails: pli9@bwh.harvard.edu and hmd.azami@gmail.com

% Example:
%   x = rand(1,1000);
%   FuzEn_MFs(x, 2, 'Bell_shaped', [0.1414*std(x) 2], 0, 1)

% Default values
if nargin < 6
    tau = 1;
end
if nargin < 5
    local = 0;
end
if nargin < 4
    rn = 0.2 * std(ts);
end

% Parse inputs
N = length(ts);

% Reconstruction
indm = hankel(1:N-m*tau, N-m*tau:N-tau);    % Indexing elements for dim-m
indm = indm(:, 1:tau:end);
ym = ts(indm);

inda = hankel(1:N-m*tau, N-m*tau:N);        % For dim-m+1
inda = inda(:, 1:tau:end);
ya = ts(inda);

% Local normalization
if local
    ym = ym - mean(ym, 2) * ones(1, m);
    ya = ya - mean(ya, 2) * ones(1, m+1);
end



% Calculate distances using 'chebyshev' metric
cheb_m = pdist(ym, 'chebychev');
cheb_a = pdist(ya, 'chebychev');

% Apply membership functions
cm = feval(mf, cheb_m, rn);
ca = feval(mf, cheb_a, rn);

% Calculate fuzzy entropy
entr = -log(sum(ca) / sum(cm));

end

% Membership functions
function c = Triangular(dist, rn)
c = zeros(size(dist));
c(dist <= rn) = 1 - dist(dist <= rn) ./ rn;
end

function c = Trapezoidal(dist, rn)
c = zeros(size(dist));
c(dist <= rn) = 1;
c(dist <= 2*rn & dist > rn) = 2 - dist(dist <= 2*rn & dist > rn) ./ rn;
end

function c = Z_shaped(dist, rn)
c = zeros(size(dist));
r1 = dist <= rn;
r2 = dist > rn & dist <= 1.5*rn;
r3 = dist > 1.5*rn & dist <= 2*rn;
c(r1) = 1;
c(r2) = 1 - 2 * ((dist(r2) - rn) ./ rn).^2;
c(r3) = 2 * ((dist(r3) - 2*rn) ./ rn).^2;
end

function c = Bell_shaped(dist, rn)
c = 1 ./ (1 + abs(dist ./ rn(1)).^(2 * rn(2)));
end

function c = Gaussian(dist, rn)
c = exp(-(dist ./ (sqrt(2) * rn)).^2);
end

function c = Constant_Gaussian(dist, rn)
c = ones(size(dist));
c(dist > rn) = exp(-log(2) .* ((dist(dist > rn) - rn) ./ rn).^2);
end

function c = Exponential(dist, rn)
c = exp(-dist.^rn(2) ./ rn(1));
end
