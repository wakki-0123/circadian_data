function [XAp, Phi] = XApEn(Sig1, Sig2, varargin)
% XApEn  estimates the cross-approximate entropy between two univariate  data sequences.
%
%   [XAp, Phi] = XApEn(Sig1, Sig2)
%
%   Returns the cross-approximate entropy estimates (``XAp``) and the log-average
%   number of matched vectors (``Phi``) for ``m`` = [0,1,2], estimated for the data
%   sequences contained in ``Sig1`` and ``Sig2`` using the default parameters:
%   embedding dimension = 2, time delay = 1,
%   radius distance threshold = 0.2*SDpooled(``Sig1``,``Sig2``),  logarithm = natural
%
%   * NOTE: ``XApEn`` is direction-dependent. Thus, the ``Sig1`` is used as
%   the template data sequence, and ``Sig2`` is the matching sequence.
%
%   [XAp, Phi] = XApEn(Sig1, Sig2, name, value, ...)
%
%   Returns the cross-approximate entropy estimates (``XAp``) between the data
%   sequences contained in ``Sig1`` and ``Sig2`` using the specified name/value pair arguments:
%
%      * ``m``     - Embedding Dimension, a positive integer   [default: 2]
%      * ``tau``   - Time Delay, a positive integer        [default: 1]
%      * ``r``     - Radius Distance Threshold, a positive scalar [default: 0.2*SDpooled(``Sig1``,``Sig2``)]
%      * ``Logx``  - Logarithm base, a positive scalar     [default: natural]
%
%   See also:
%       XSampEn, XFuzzEn, XMSEn, ApEn, SampEn, MSEn
%
%   References:
%     [1] Steven Pincus and Burton H. Singer,
%           "Randomness and degrees of irregularity."
%           Proceedings of the National Academy of Sciences
%           93.5 (1996): 2083-2088.
%
%     [2] Steven Pincus,
%           "Assessing serial irregularity and its implications for health."
%           Annals of the New York Academy of Sciences
%           954.1 (2001): 245-267.
%

narginchk(2,10)
p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk1 = @(x) isnumeric(x) && length(x)>=10 && isvector(x);
% Chk11 = @(x) (isnumeric(x) && isvector(x) && numel(x)>=10) || isempty(x);
Chk2 = @(x) isscalar(x) && (x > 0);
addRequired(p,'Sig1',Chk1);
addRequired(p,'Sig2',Chk1);
S1 = Sig1(:); S2 = Sig2(:);
N1 = length(S1); N2 = length(S2);

addParameter(p,'m',2,Chk);
addParameter(p,'tau',1,Chk);
addParameter(p,'r',0.2*sqrt((var(S1,1)*(N1-1) + var(S2,1)*(N2-1))/(N1+N2-1)),Chk2);
addParameter(p,'Logx',exp(1),Chk2);
parse(p, Sig1, Sig2, varargin{:})
m = p.Results.m; tau = p.Results.tau;
r = p.Results.r; Logx = p.Results.Logx;

Counter = cast(1*(abs(S1 - S2') <= r),"int8");
M = [m*ones(1,N1-(m*tau)) repelem((m-1):-1:1,tau)];
XAp = zeros(1,m); Phi = XAp;
for n = 1:N1-tau
    ix = find(Counter(n,:)==1);
    for k = 1:M(n)
        ix(ix + (k*tau) > N2) = [];
        if isempty(ix)
            break
        end
        p1 = repmat(S1(n:tau:n+k*tau)',length(ix),1);
        p2 = S2(ix+(0:tau:k*tau)')';
        ix = ix(any(max(abs(p1 - p2),[],2) <= r,2));
        Counter(n, ix) = Counter(n, ix) + 1;
    end
end

Phi(1) = 0; %(log(N1)/log(Logx))/N1;
Temp = sum(Counter>0)/N1; Temp(Temp==0) = [];
Phi(2) = mean(log(Temp)/log(Logx));
XAp(1) = Phi(1) - Phi(2);
for k = 0:m-1
    ai = sum(Counter>k+1)/(N1-(k+1)*tau);
    bi = sum(Counter>k)/(N1-(k*tau));
    ai(ai==0) = [];
    bi(bi==0) = [];
    Phi(k+3) = sum(log(ai)/log(Logx))/(N1-(k+1)*tau);
    XAp(k+2)= sum(log(bi)/log(Logx))/(N1-(k*tau))- Phi(k+3);
end
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub