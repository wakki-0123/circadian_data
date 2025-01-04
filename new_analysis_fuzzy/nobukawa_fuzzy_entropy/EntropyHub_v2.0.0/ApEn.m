function [Ap, Phi] = ApEn(Sig, varargin)
% ApEn  estimates the approximate entropy of a univariate data sequence.
%
%   [Ap, Phi] = ApEn(Sig) 
%
%   Returns the approximate entropy estimates (``Ap``) and the log-average number of 
%   matched vectors (``Phi``) for ``m`` = [0, 1, 2], estimated from the data 
%   sequence (``Sig``) using the default parameters:
%   embedding dimension = 2, time delay = 1,
%   radius distance threshold = 0.2*SD(``Sig``), logarithm = natural
%
%   [Ap, Phi] = ApEn(Sig, name, value, ...)
%
%   Returns the approximate entropy estimates (``Ap``) of the data sequence (``Sig``)
%   for dimensions = [0, 1, ..., ``m``] using the specified name/value pair arguments:
% 
%       * ``m``     - Embedding Dimension, a positive integer
%       * ``tau``   - Time Delay, a positive integer
%       * ``r``     - Radius Distance Threshold, a positive scalar  
%       * ``Logx``  - Logarithm base, a positive scalar  
%
%   See also: 
%       XApEn, SampEn, MSEn, FuzzEn, PermEn, CondEn, DispEn
%   
%   References:
%     [1] Steven M. Pincus, 
%           "Approximate entropy as a measure of system complexity." 
%           Proceedings of the National Academy of Sciences 
%           88.6 (1991): 2297-2301.
% 

narginchk(1,9)
if size(Sig,1) == 1
    Sig = Sig';
end
Sig = squeeze(Sig);

p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isscalar(x) && (x > 0);
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'m',2,Chk);
addParameter(p,'tau',1,Chk);
addParameter(p,'r',.2*std(Sig,1),Chk2);
addParameter(p,'Logx',exp(1),Chk2);
parse(p,Sig,varargin{:})
m = p.Results.m; tau = p.Results.tau; 
r = p.Results.r; Logx = p.Results.Logx; 

N = length(Sig);
Counter = 1*(abs(Sig - Sig') <= r);
M = [m*ones(1,N-(m*tau)) repelem((m-1):-1:1,tau)];
Ap = zeros(1,m); Phi = Ap;
for n = 1:N-tau
    ix = find(Counter(n,:)==1);
    for k = 1:M(n)
        ix(ix + (k*tau) > N) = [];  
        p1 = repmat(Sig(n:tau:n+k*tau)',length(ix),1);
        p2 = Sig(ix+(0:tau:k*tau)')';
        ix = ix(max(abs(p1 - p2),[],2) <= r);
        Counter(n, ix) = Counter(n, ix) + 1;   
    end
end

Phi(1) = 0; %(log(N)/log(Logx))/N; 
Phi(2) = mean(log(sum(Counter>0)/N)/log(Logx));
Ap(1) = Phi(1) - Phi(2);
for k = 0:m-1
    ai = sum(Counter>k+1)/(N-(k+1)*tau);
    bi = sum(Counter>k)/(N-(k*tau));
    ai(ai==0) = [];
    bi(bi==0) = [];  
    Phi(k+3) = sum(log(ai)/log(Logx))/(N-(k+1)*tau);
    Ap(k+2)= sum(log(bi)/log(Logx))/(N-(k*tau))- Phi(k+3); 
end

end
%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub