function [Div, CDEn, Bm] = DivEn(Sig, varargin) 
% DivEn  estimates the diversity entropy of a univariate data sequence.
%
%   [Div, CDEn, Bm] = DivEn(Sig)
%
%   Returns the diversity entropy (``Div``), the cumulative diversity entropy (``CDEn``),
%   and the corresponding probabilities (``Bm``) estimated from the data sequence (``Sig``)
%   using the default parameters:   embedding dimension = 2, time delay = 1,
%   # bins = 5,  logarithm = natural,
%
%   [Div, CDEn, Bm] = DivEn(Sig, name, value, ...)
%
%   Returns the diversity entropy (``Div``) estimated from the data
%   sequence (``Sig``) using the specified name/value pair arguments:
%       * ``m``     - Embedding Dimension, an integer > 1
%       * ``tau``   - Time Delay, a positive integer
%       * ``r``     - Histogram bins #: either
%                        * an integer [``r`` > 1] representing the number of bins
%                        * a vector array of 3 or more increasing values in range [-1 1] representing the bin edges including the rightmost edge.
%       * ``Logx``  - Logarithm base, a positive scalar (enter 0 for natural log)
%
%   See also:
%      CoSiEn, PhasEn, SlopEn, GridEn, MSEn
%
%    References:
%       [1] X. Wang, S. Si and Y. Li,
%           "Multiscale Diversity Entropy: A Novel Dynamical Measure for Fault
%           Diagnosis of Rotating Machinery,"
%           IEEE Transactions on Industrial Informatics,
%           vol. 17, no. 8, pp. 5419-5429, Aug. 2021
%
%       [2] Y. Wang, M. Liu, Y. Guo, F. Shu, C. Chen and W. Chen,
%           "Cumulative Diversity Pattern Entropy (CDEn): A High-Performance,
%           Almost-Parameter-Free Complexity Estimator for Nonstationary Time Series,"
%           IEEE Transactions on Industrial Informatics
%           vol. 19, no. 9, pp. 9642-9653, Sept. 2023
%

narginchk(1,9)
Sig = squeeze(Sig);
N = numel(Sig);
p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && ((isscalar(x) && x>1) ||  ...
    (isvector(x) && length(x)>2 && min(x)>=-1 && max(x)<=1 && min(diff(x))>0));
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'m',2,Chk);
addParameter(p,'tau',1,Chk);
addParameter(p,'r',5,Chk2);
addParameter(p,'Logx',exp(1),@(x) isscalar(x) && (x > 0));
parse(p,Sig,varargin{:})
m = p.Results.m; tau = p.Results.tau; 
r = p.Results.r; Logx = p.Results.Logx; 
if Logx == 0
    Logx = exp(1);
end
if isscalar(r)
    r = linspace(-1,1,r+1);
end

Nx = N - (m-1)*tau;
Zm = zeros(Nx,m);
for n = 1:m
    Zm(:,n) = Sig((n-1)*tau+1:Nx+(n-1)*tau);
end

Num = sum(Zm(1:end-1,:).*Zm(2:end,:),2)';
Den = vecnorm(Zm(2:end,:)',2).*vecnorm(Zm(1:end-1,:)',2);
Di = Num./Den;
Bm = histcounts(Di, r); %/(Nx-1);
Bm = Bm(Bm>0)/sum(Bm);

if round(sum(Bm),6) ~= 1
    warning("Potential error is probability estimation!\n" + ...
        "Sum(Pi) == %d", round(sum(Bm),6))
end
r = length(r)-1;

Pj = 1 - cumsum(Bm);     
Pj = (Pj/sum(Pj));
Pj = Pj(1:end-1);
CDEn = -sum(Pj.*log(Pj)/log(Logx))/(log(r)/log(Logx));
Div = -sum(Bm.*log(Bm)/log(Logx))/(log(r)/log(Logx));

end


%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub