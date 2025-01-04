function [Incr] = IncrEn(Sig, varargin)
% IncrEn  estimates the increment entropy of a univariate data sequence.
%
%   [Incr] = IncrEn(Sig) 
% 
%   Returns the increment entropy (``Incr``) estimated from the data sequence 
%   (``Sig``) using the default parameters: 
%   embedding dimension = 2, time delay = 1, quantifying resolution = 4,
%   logarithm = base 2,
%
%   [Incr] = IncrEn(Sig, name, value, ...)
% 
%   Returns the increment entropy (``Incr``) estimated from the data sequence 
%   (``Sig``) using the specified name/value pair arguments:
% 
%      * ``m``     - Embedding Dimension, an integer > 1
%      * ``tau``   - Time Delay, a positive integer
%      * ``R``     - Quantifying resolution, a positive integer
%      * ``Logx``  - Logarithm base, a positive scalar (enter 0 for natural log) 
%      * ``Norm``  - Normalisation of IncrEn value, a boolean:
%                * [false]  no normalisation - default
%                * [true]   normalises w.r.t embedding dimension (m-1). 
%
%   See also:
%       PermEn, SyDyEn, MSEn
%   
%   References:
%      [1] Xiaofeng Liu, et al.,
%           "Increment entropy as a measure of complexity for time series."
%           Entropy
%           18.1 (2016): 22.1.
% 
%      ***   "Correction on Liu, X.; Jiang, A.; Xu, N.; Xue, J. - Increment 
%           Entropy as a Measure of Complexity for Time Series,
%           Entropy 2016, 18, 22." 
%           Entropy 
%           18.4 (2016): 133.
% 
%      [2] Xiaofeng Liu, et al.,
%           "Appropriate use of the increment entropy for 
%           electrophysiological time series." 
%           Computers in biology and medicine 
%           95 (2018): 13-23.
% 


narginchk(1,11)

p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
Chkx = @(x) isnumeric(x) && isscalar(x) && (x > 1) && (mod(x,1)==0);
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'m',2,Chkx);
addParameter(p,'tau',1,Chk);
addParameter(p,'R',4,@(x) isnumeric(x) && (x > 0) && (mod(x,1)==0));
addParameter(p,'Logx',2,Chk2);
addParameter(p,'Norm',false,@(x) islogical(x));
parse(p,Sig,varargin{:})
m = p.Results.m; tau = p.Results.tau; R = p.Results.R;
Logx = p.Results.Logx; Norm = p.Results.Norm;

if Logx == 0
    Logx = exp(1);
end
Vi = diff(squeeze(Sig));
N = length(Vi)-((m-1)*tau);
Vk = zeros(N,m);
for k = 1:m
    Vk(:,k) = Vi(1+(k-1)*tau:N+(k-1)*tau);
end

Sk = sign(Vk);
Temp = std(Vk,[],2);
Qk = min(R, floor((abs(Vk)*R)./repmat(Temp,1,m)));
Qk(any(Temp==0,2),:) = 0;
Wk = Sk.*Qk;
Px = unique(Wk,'rows');
Counter = zeros(1,size(Px,1));
for k = 1:size(Px,1) 
    Counter(k) = sum(~any(Wk - Px(k,:),2));
end
Ppi = Counter/N;

if size(Px,1) > (2*R + 1)^m
    warning('Error with probability estimation')
elseif round(sum(Ppi),3) ~= 1
    warning('Error with probability estimation')
end
Incr = -sum(Ppi.*(log(Ppi)/log(Logx)));
if Norm
    Incr = Incr/(m-1);
end
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub