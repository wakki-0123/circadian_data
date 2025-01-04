function [XK2, Ci] = XK2En(Sig1, Sig2, varargin) 
% XK2En  estimates the cross-Kolmogorov (K2) entropy between two univariate data sequences.
%
%   [XK2, Ci] = XK2En(Sig1, Sig2)
% 
%   Returns the cross-Kolmogorov entropy estimates (``XK2``) and the correlation
%   integrals (``Ci``) for ``m`` = [1,2] estimated between the data sequences 
%   contained in ``Sig1`` and ``Sig2`` using the default parameters: 
%   embedding dimension = 2, time delay = 1, 
%   distance threshold (``r``) = 0.2*SDpooled(``Sig1``,``Sig2``), logarithm = natural
%
%   [XK2, Ci] = XK2En(Sig, name, value, ...)
% 
%   Returns the cross-Kolmogorov entropy estimates (``XK2``) estimated between
%   the data sequences contained in ``Sig1`` and ``Sig2`` using the specified name/value
%   pair arguments:
% 
%       * ``m``     - Embedding Dimension, a positive integer [default: 2]
%       * ``tau``   - Time Delay, a positive integer         [default: 1]
%       * ``r``     - Radius Distance Threshold, a positive scalar [default: 0.2*SDpooled(``Sig1``,``Sig2``)]
%       * ``Logx``  - Logarithm base, a positive scalar      [default: natural]
%
%   See also:
%      XSampEn, XFuzzEn, XApEn, K2En, XMSEn, XDistEn.
% 
%   References:
%     [1]   Matthew W. Flood,
%               "XK2En - EntropyHub Project"
%               (2021) https://github.com/MattWillFlood/EntropyHub
% 

narginchk(2,10)
p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isscalar(x) && (x > 0) && isnumeric(x);
Chk1 = @(x) isnumeric(x) && isvector(x) && numel(x)>=10;
addRequired(p,'Sig1',Chk1);
addRequired(p,'Sig2',Chk1);
S1 = Sig1(:);  S2 = Sig2(:);
N1 = length(S1);  N2 = length(S2);

addParameter(p,'m',2,Chk);
addParameter(p,'tau',1,Chk);
addParameter(p,'r',0.2*sqrt((var(S1,1)*(N1-1) + var(S2,1)*(N2-1))/(N1+N2-1)),Chk2);
addParameter(p,'Logx',exp(1),Chk2);
parse(p,Sig1, Sig2, varargin{:})
m = p.Results.m; tau = p.Results.tau; 
r = p.Results.r; Logx = p.Results.Logx; 

m   = m+1;
Zm1 = zeros(N1,m);
Zm2 = zeros(N2,m); 
Ci = zeros(1,m);
for n = 1:m
    Nx = N1-(n-1)*tau;
    Zm1(1:Nx,n) = S1((n-1)*tau + 1:N1); 
    Ny = N2-(n-1)*tau;
    Zm2(1:Ny,n) = S2((n-1)*tau + 1:N2);   
    Norm = zeros(Nx,Ny);    
    for k = 1:Nx
        Temp = repmat(Zm1(k,1:n),Ny,1) - Zm2(1:Ny,1:n);
        Norm(k,:) = sqrt(sum(Temp.*Temp,2)); 
    end
    Ci(n) = mean(Norm(:) < r);     
end
 
XK2 = (log(Ci(1:m-1)./Ci(2:m))/log(Logx))/tau;
XK2(isinf(XK2)) = NaN;
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub