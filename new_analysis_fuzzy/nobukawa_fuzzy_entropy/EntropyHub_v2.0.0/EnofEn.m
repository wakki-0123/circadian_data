function [EoE, AvEn, S2] = EnofEn(Sig, varargin) 
% EnofEn  estimates the entropy of entropy from a univariate data sequence.
%
%   [EoE, AvEn, S2] = EnofEn(Sig) 
% 
%   Returns the entropy of entropy (``EoE``), the average Shannon entropy 
%   (``AvEn``), and the number of levels (``S2``) across all windows 
%   estimated from the data sequence (``Sig``) using the default parameters: 
%   window length (samples) = 10, slices = 10, logarithm = natural
%   heartbeat interval range (xmin, xmax) = [min(Sig) max(Sig)]
%
%   [EoE, AvEn, S2] = EnofEn(Sig, name, value, ...)
% 
%   Returns the entropy of entropy (``EoE``) estimated from the data sequence (``Sig``)  
%   using the specified name/value pair arguments:
% 
%      * ``tau``    - Window length, an integer > 1
%      * ``S``      - Number of slices, an integer > 1
%      * ``Xrange`` - The min and max heartbeat interval,
%                     a two-element vector where X(1) < X(2)
%      * ``Logx``   - Logarithm base, a positive scalar  
% 
%   See also:
%       SampEn, MSEn
%   
%   References:
%      [1] Chang Francis Hsu, et al.,
%           "Entropy of entropy: Measurement of dynamical complexity for
%           biological systems." 
%           Entropy 
%           19.10 (2017): 550.
% 


narginchk(1,9)
Sig = squeeze(Sig);

p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x)  && (x<length(Sig)) && (x>1) && (mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (x > 0);
Chk3 = @(x) isnumeric(x) && (x > 1) && (mod(x,1)==0);
Chk4 = @(x) isnumeric(x) && (numel(x)==2) && (diff(x)>=0);

addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'S',10,Chk3);
addParameter(p,'tau',10,Chk);
addParameter(p,'Xrange',[min(Sig) max(Sig)],Chk4);
addParameter(p,'Logx',exp(1),Chk2);
parse(p,Sig,varargin{:})
tau = p.Results.tau; Logx = p.Results.Logx; S = p.Results.S; 
Xrange = p.Results.Xrange;

Wn = floor(length(Sig)/tau);
Wj = reshape(Sig(1:Wn*tau),tau,Wn)';
Yj = zeros(1,Wn);
Edges = linspace(Xrange(1),Xrange(2),S+1); % Edges = linspace(min(Sig),max(Sig),S(1)+1);

for n = 1:Wn
    Temp = histcounts(Wj(n,:),Edges)/tau;
    Temp(Temp==0) = [];
    Yj(n) = -sum(Temp.*(log(Temp)/log(Logx)));
end

AvEn = sum(Yj)/Wn;
% Edges = linspace(min(Yj),max(Yj),S+1);
% Pjl = histcounts(Yj,Edges)/Wn;
% Pjl(Pjl==0) = [];
[~,~,Tempy] = unique(round(Yj,12));
Pjl = accumarray(Tempy,1)/Wn;
if round(sum(Pjl),5) ~= 1
    warning('Possible error estimating probabilities')
end
S2 = length(Pjl);
EoE = -sum(Pjl.*(log(Pjl)/log(Logx)));
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub