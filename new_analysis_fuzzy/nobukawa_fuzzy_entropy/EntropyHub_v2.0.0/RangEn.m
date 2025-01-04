function [Rangx, A, B] =  RangEn(Sig, varargin)
%  RangEn  estimates the range entropy of a univariate data sequence.
%
%   [Rangx, A, B] = RangEn(Sig)
%
%   Returns the range entropy estimate (``Rangx``) and the number of matched state
%   vectors (``m: B``, ``m+1: A``) estimated from the data sequence (``Sig``)
%   using the sample entropy algorithm and the following default parameters:
%   embedding dimension = 2, time delay = 1, radius threshold = 0.2, logarithm = natural.
%
%   Rangx, A, B = RangEn(Sig, name, value, ...)
%
%   Returns the range entropy estimates (``Rangx``) for dimensions = ``m``
%   estimated for the data sequence (``Sig``) using the specified name-value arguments:
%       * ``m``         - Embedding Dimension, a positive integer
%       * ``tau``       - Time Delay, a positive integer
%       * ``r``         - Radius Distance Threshold, a positive scalar between 0 and 1
%       * ``Methodx``   - Base entropy method, either 'SampEn' [default] or 'ApEn'
%       * ``Logx``      - Logarithm base, a positive scalar
%
%   See also:
%       ``ApEn``, ``SampEn``, ``FuzzEn``,  ``MSEn``
%
%   References:
%    [1] Omidvarnia, Amir, et al.
%        "Range entropy: A bridge between signal complexity and self-similarity"
%        Entropy
%        20.12 (2018): 962.
%
%    [2] Joshua S Richman and J. Randall Moorman.
%        "Physiological time-series analysis using approximate entropy
%        and sample entropy."
%        American Journal of Physiology-Heart and Circulatory Physiology
%        2000
%

narginchk(1,11)
Sig = squeeze(Sig);
if size(Sig,1) == 1
    Sig = Sig';
end
N = size(Sig,1);
p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
addRequired(p,'Sig', @(x) isnumeric(x) && isvector(x) && (length(x)>=10));
addParameter(p,'m',2,Chk);
addParameter(p,'tau',1,Chk);
addParameter(p,'r', 0.2, @(x) isnumeric(x) && isscalar(x) && x>=0 && x<=1);
addParameter(p,'Methodx','sampen',@(x) any(validatestring(lower(x),{'sampen','apen'})));
addParameter(p,'Logx',exp(1),@(x) isscalar(x) && (x > 0));
parse(p,Sig,varargin{:})
m = p.Results.m; tau = p.Results.tau;
r = p.Results.r; Logx = p.Results.Logx;
Methodx = p.Results.Methodx;
if Logx == 0
    Logx = exp(1);
end

switch lower(Methodx)
    case 'sampen'
        Nx = N - m*tau;
        Sx = zeros(Nx,m+1);
        for k = 1:m+1
            Sx(:,k) = Sig(1+(k-1)*tau:Nx + (k-1)*tau);
        end
        A = zeros(Nx,1);
        B = zeros(Nx,1);
        for k = 1:(Nx - 1)
            Dxy = abs(repmat(Sx(k,1:end-1),Nx-k,1) - Sx(k+1:end,1:end-1));
            Mx = max(Dxy,[],2);
            Mn = min(Dxy,[],2);
            RR = (Mx - Mn)./(Mx + Mn) <= r;
            B(k) = sum(RR);

            if B(k)>0
                Temp = Sx(k+1:end,:);
                Dxy2 = abs(repmat(Sx(k,:),B(k),1) - Temp(RR,:));
                Mx = max(Dxy2,[],2);
                Mn = min(Dxy2,[],2);
                RR2 = (Mx - Mn)./(Mx + Mn) <= r;
                A(k) = sum(RR2);
            end
        end
        Rangx = -log(sum(A)/sum(B))/log(Logx);

    case 'apen'
        Nx = N - (m-1)*tau;
        Sx = zeros(Nx,m);
        for k = 1:m
            Sx(:,k) = Sig((k-1)*tau + 1:Nx + (k-1)*tau);
        end
        Sx = [Sx [Sig(m*tau + 1:end); zeros(tau,1)]];

        B = zeros(Nx,1);
        A = zeros(Nx-tau,1);
        for k = 1:Nx
            %Dxy = abs(repmat(Sx(k,1:end-1),Nx,1) - Sx(:,1:end-1));
            Dxy = abs(repmat(Sx(k,1:m),Nx,1) - Sx(:,1:m));
            Mx = max(Dxy,[],2);
            Mn = min(Dxy,[],2);
            RR = (Mx - Mn)./(Mx + Mn) <= r;
            B(k) = sum(RR);        
            if k <= (Nx - tau)
                RR(end-tau+1:end) = false;
                Dxy2 = abs(repmat(Sx(k,:),sum(RR),1) - Sx(RR,:));
                Mx2 = max(Dxy2,[],2);
                Mn2 = min(Dxy2,[],2);
                RR2 = (Mx2 - Mn2)./(Mx2 + Mn2) <= r;
                A(k) = sum(RR2);
            end
        end
        Ax = mean(log(A/(Nx-tau))/log(Logx));
        Bx = mean(log(B/Nx)/log(Logx));
        Rangx = Bx - Ax;

    otherwise
        error("Methodx must be either 'ApEn' or 'SampEn'")
end
end

% Copyright 2024 Matthew W. Flood, EntropyHub
% For Terms of Use see https://github.com/MattWillFlood/EntropyHub