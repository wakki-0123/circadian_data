function [Dist, Ppi] = DistEn(Sig, varargin)
% DistEn  estimates the distribution entropy of a univariate data sequence.
%
%   [Dist, Ppi] = DistEn(Sig) 
% 
%   Returns the distribution entropy estimate (``Dist``) and the
%   corresponding distribution probabilities (``Ppi``) estimated from the data 
%   sequence (``Sig``) using the default  parameters: 
%   embedding dimension = 2, time delay = 1, binning method = ``'Sturges'``,
%   logarithm = base 2, normalisation = w.r.t # of histogram bins
%
%   [Dist, Ppi] = DistEn(Sig, name, value, ...)
% 
%   Returns the distribution entropy estimate (``Dist``) estimated from the data
%   sequence (``Sig``) using the specified name/value pair arguments:
% 
%       * ``m``     - Embedding Dimension, a positive integer
%       * ``tau``   - Time Delay, a positive integer
%       * ``Bins``  - Histogram bin selection method for distance distribution, either
%         an integer > 1 indicating the number of bins, or one of the following strings 
%         {``'sturges'``, ``'sqrt'``, ``'rice'``, ``'doanes'``} [default: ``'sturges'``]
%       * ``Logx``  - Logarithm base, a positive scalar (enter 0 for natural log) 
%       * ``Norm``  - Normalisation of ``Dist`` value, a boolean:
%                 - [false]  no normalisation.
%                 - [true]   normalises w.r.t # of histogram bins (default)
% 
%   See also:
%       XDistEn, DistEn2D, MSEn, K2En
%   
%   References:
%     [1] Li, Peng, et al.,
%           "Assessing the complexity of short-term heartbeat interval 
%           series by distribution entropy." 
%           Medical & biological engineering & computing 
%           53.1 (2015): 77-87. 
% 


narginchk(1,11)
Sig = squeeze(Sig);

p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isscalar(x) && (x >= 0);
Chk3 = {'sturges','sqrt','rice','doanes'};
Chk4 = @(x) (isnumeric(x) && (length(x)==1) && (mod(x,1)==0)) ...
    || (ischar(x) && any(validatestring(lower(x),Chk3)));
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'m',2,Chk);
addParameter(p,'tau',1,Chk);
addParameter(p,'Bins','sturges',Chk4);
addParameter(p,'Logx',2,Chk2);
addParameter(p,'Norm',true,@(x) islogical(x));
parse(p,Sig,varargin{:})
m = p.Results.m; tau = p.Results.tau; Bins = p.Results.Bins; 
Logx = p.Results.Logx; Norm = p.Results.Norm;

if Logx == 0
    Logx = exp(1);
end
Nx = length(Sig) - ((m-1)*tau);
Zm = zeros(Nx,m);
for n = 1:m
    Zm(:,n) = Sig((n-1)*tau + 1:Nx+(n-1)*tau);
end

DistMat = zeros(1,Nx*(Nx-1)/2);
for k = 1:Nx-1
    Ix = [((k-1)*(Nx - k/2)+1), k*(Nx-((k+1)/2))];
    DistMat(Ix(1):Ix(2)) = max(abs(repmat(Zm(k,:),Nx-k,1) - Zm(k+1:end,:)),[],2);
end

Ny = length(DistMat);
if ischar(Bins)
    switch lower(Bins)
        case 'sturges'
            Bx = ceil(log2(Ny) + 1);
        case 'rice'
            Bx = ceil(2*(Ny^(1/3)));
        case 'sqrt'
            Bx = ceil(sqrt(Ny));
        case 'doanes'
            sigma = sqrt(6*(Ny-2)/((Ny+1)*(Ny+3)));
            Bx = ceil(1+log2(Ny)+log2(1+abs(skewness(DistMat)/sigma)));
        otherwise 
            error('Please enter a valid binning method')
    end
else    
    Bx = Bins;
end
By = linspace(min(DistMat),max(DistMat),Bx+1);
Ppi = histcounts(DistMat,By)/Ny;

if round(sum(Ppi),6) ~= 1
    warning('Potential error estimating probabilities (p=%d).', sum(Ppi))
    Ppi(Ppi==0)=[];
elseif any(Ppi==0)
    fprintf('Note: %d/%d bins were empty',sum(Ppi==0),numel(Ppi));
    Ppi(Ppi==0)=[];
end
Dist = -sum(Ppi.*log(Ppi)/log(Logx));
if Norm == 1
    Dist = Dist/(log(Bx)/log(Logx));
end
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub