function [MCoSi, Bm] = MvCoSiEn(Data, varargin)
% MvCoSiEn  estimates the multivariate cosine similarity entropy of a multivariate dataset.
%
%   [MCoSi, Bm] = MvCoSiEn(Data) 
% 
%   Returns the multivariate cosine similarity entropy estimate (``MCoSi``)
%   and the corresponding global probabilities (``Bm``) estimated for the 
%   M multivariate sequences in ``Data`` using the default parameters: 
%   embedding dimension = 2*ones(M,1), time delay = ones(M,1), angular threshold = 0.1,
%   logarithm = 2, data normalization = none, 
% 
%   .. note::
%       To maximize the number of points in the embedding process, this algorithm 
%       uses N-max(m*tau) delay vectors and **not** N-max(m)*max(tau) as employed 
%       in [1][2][3].
% 
%   [MCoSi, Bm] = MvCoSiEn(Data, name, value, ...)
% 
%   Returns the multivariate cosine similarity entropy estimates (``MCoSi``) estimated
%   from the M multivariate data sequences in ``Data`` using the specified name/value pair arguments:
% 
%       :Data:  - Multivariate dataset, NxM matrix of N (>10) observations (rows) and M (cols) univariate data sequences 
%       :m:     - Embedding Dimension, a vector of M positive integers
%       :tau:   - Time Delay, a vector of M positive integers
%       :r:     - Angular threshold, a value in range [0 < ``r`` < 1]
%       :Logx:  - Logarithm base, a positive scalar (enter 0 for natural log) 
%       :Norm:  - Normalisation of ``Data``, one of the following integers
%                     :0:  no normalisation (default)
%                     :1:  remove median(``Data``) to get zero-median series
%                     :2:  remove mean(``Data``) to get zero-mean series
%                     :3:  normalises each sequence in ``Data`` to unit variance and zero mean
%                     :4:  normalises each sequence in ``Data`` values to range [-1 1]
%
%   See also:
%       CoSiEn, MvMSEn, MvSampEn, MvDispEn, MvFuzzEn, MvPermEn, cMvMSEn
%   
%   References:
%      [1] H. Xiao, T. Chanwimalueang and D. P. Mandic, 
%           "Multivariate Multiscale Cosine Similarity Entropy" 
%           IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP),
%           pp. 5997-6001, doi: 10.1109/ICASSP43922.2022.9747282.
% 
%      [2] Xiao, H.; Chanwimalueang, T.; Mandic, D.P., 
%           "Multivariate Multiscale Cosine Similarity Entropy and Its 
%           Application to Examine Circularity Properties in Division Algebras."
%           Entropy 2022, 24, 1287. 
% 
%      [3] Ahmed Mosabber Uddin, Danilo P. Mandic
%           "Multivariate multiscale entropy: A tool for complexity
%           analysis of multichannel data."
%           Physical Review E 84.6 (2011): 061918.
% 
%      [4] Theerasak Chanwimalueang and Danilo Mandic,
%           "Cosine similarity entropy: Self-correlation-based complexity
%           analysis of dynamical systems."
%           Entropy 
%           19.12 (2017): 652.
% 

narginchk(1,11)
p = inputParser;
Data = squeeze(Data);
Dn = size(Data,2);
N = size(Data,1);
Chk = @(x) isnumeric(x) && isvector(x) && (length(x)==Dn) && (min(x)>0) && all(mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (x >= 0);

addRequired(p,'Data',@(x) isnumeric(x) && ismatrix(x) && (Dn>1) && N>10);
addParameter(p,'m',2*ones(Dn,1),Chk);
addParameter(p,'tau',ones(Dn,1),Chk);
addParameter(p,'r',.1,@(x) isnumeric(x) && isscalar(x) && (x > 0) && (x < 1));
addParameter(p,'Logx',2,Chk2);
addParameter(p,'Norm',0,@(x) mod(x,1)==0 && ismember(x,[0:4]));
parse(p,Data,varargin{:})
m = p.Results.m(:); tau = p.Results.tau(:); 
r = p.Results.r(:); Logx = p.Results.Logx; 
Norm = p.Results.Norm;

if Logx == 0;     Logx = exp(1); end
if Norm == 1
    Xi = Data - median(Data,1);
elseif Norm == 2
    Xi = Data - mean(Data,1);
elseif Norm == 3
    Xi = (Data - mean(Data,1))./std(Data,1,1);
elseif Norm == 4
    Xi = (2*(Data - min(Data,[],1))./range(Data,1)) - 1;
else
    Xi = Data;
end

Nx = N - max((m-1).*tau);
Zm = zeros(Nx,sum(m));
q = 1;
for k = 1:Dn
    for p=1:m(k)
        Zm(:,q) = Xi(1+(p-1)*tau(k):Nx+(p-1)*tau(k),  k);
        q = q+ 1;
    end    
end

Num = Zm*Zm'; 
Mag = sqrt(diag(Num));
Den = Mag*Mag';
AngDis = acos(Num./Den)/pi;
if max(imag(AngDis(:))) < (10^-5) %max(max(imag(AngDis))) < (10^-5)
    Bm = sum(sum(triu(round(AngDis,6) < r,1)))/(Nx*(Nx-1)/2);
else
    Bm = sum(sum(triu(real(AngDis) < r,1)))/(Nx*(Nx-1)/2);
    warning('Complex values ignored')
end
if Bm == 1 || Bm == 0
    MCoSi = nan;
else
    MCoSi = -(Bm*log(Bm)/log(Logx)) - ((1-Bm)*log(1-Bm)/log(Logx));
end

end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub