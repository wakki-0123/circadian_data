function [CoSi, Bm] = CoSiEn(Sig, varargin)
% CoSiEn  estimates the cosine similarity entropy of a univariate data sequence.
%
%   [CoSi, Bm] = CoSiEn(Sig) 
% 
%   Returns the cosine similarity entropy (``CoSi``) and the corresponding
%   global probabilities (``Bm``) estimated from the data sequence (``Sig``) 
%   using the default parameters: 
%   embedding dimension = 2, time delay = 1, angular threshold = .1,
%   logarithm = base 2,
%
%   [CoSi, Bm] = CoSiEn(Sig, name, value, ...)
% 
%   Returns the cosine similarity entropy (``CoSi``) estimated from the data
%   sequence (``Sig``) using the specified name/value pair arguments:
% 
%       * ``m``     - Embedding Dimension, an integer > 1
%       * ``tau``   - Time Delay, a positive integer
%       * ``r``     - Angular threshold, a value in range [0 < ``r`` < 1]
%       * ``Logx``  - Logarithm base, a positive scalar (enter 0 for natural log) 
%       * ``Norm``  - Normalisation of ``Sig``, one of the following integers:
%               *  [0]  no normalisation - default
%               *  [1]  remove median(``Sig``) to get zero-median series
%               *  [2]  remove mean(``Sig``) to get zero-mean series
%               *  [3]  normalises ``Sig`` w.r.t. SD(``Sig``)
%               *  [4]  normalises ``Sig`` values to range [-1 1]
%
%   See also:
%       PhasEn, SlopEn, GridEn, MSEn, hMSEn
%   
%   References:
%     [1] Theerasak Chanwimalueang and Danilo Mandic,
%           "Cosine similarity entropy: Self-correlation-based complexity
%           analysis of dynamical systems."
%           Entropy 
%           19.12 (2017): 652.
% 

narginchk(1,11)
Sig = squeeze(Sig);

p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
Chk3 = @(x) isnumeric(x) && isscalar(x) && (x > 1) && (mod(x,1)==0);
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'m',2,Chk3);
addParameter(p,'tau',1,Chk);
addParameter(p,'r',.1,@(x) isnumeric(x) && (x > 0) && (x < 1));
addParameter(p,'Logx',2,Chk2);
addParameter(p,'Norm',0,@(x) ismember(x,[0:4]));
parse(p,Sig,varargin{:})
m = p.Results.m; tau = p.Results.tau; r = p.Results.r;
Logx = p.Results.Logx; Norm = p.Results.Norm;

N = length(Sig);
if Logx == 0
    Logx = exp(1);
end
if Norm == 1
    Xi = Sig - median(Sig);
elseif Norm == 2
    Xi = Sig - mean(Sig);
elseif Norm == 3
    Xi = (Sig - mean(Sig))/std(Sig,1);
elseif Norm == 4
    Xi = (2*(Sig - min(Sig))/range(Sig)) - 1;
else
    Xi = Sig;
end
Nx = N-((m-1)*tau);
Zm = zeros(Nx,m);
for n = 1:m
    Zm(:,n) = Xi((n-1)*tau+1:Nx+(n-1)*tau);
end

Num = Zm*Zm'; 
Mag = sqrt(diag(Num));
Den = Mag*Mag';
AngDis = acos(Num./Den)/pi;
if max(imag(AngDis(:))) < (10^-5) 
    Bm = sum(sum(triu(round(AngDis,6) < r,1)))/(Nx*(Nx-1)/2);
else
    Bm = sum(sum(triu(real(AngDis) < r,1)))/(Nx*(Nx-1)/2);
    warning('Complex values ignored')
end
if Bm == 1 || Bm == 0
    CoSi = nan;
else
    CoSi = -(Bm*log(Bm)/log(Logx)) - ((1-Bm)*log(1-Bm)/log(Logx));
end
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub