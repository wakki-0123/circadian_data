function [Cond, SEw, SEz] = CondEn(Sig, varargin)
% CondEn  estimates the corrected conditional entropy of a univariate data sequence.
%
%   [Cond, SEw, SEz] = CondEn(Sig) 
% 
%   Returns the corrected conditional entropy estimates (``Cond``) and the
%   corresponding Shannon entropies (``m: SEw``, ``m+1: SEz``) for ``m`` = [1,2] 
%   estimated from the data sequence (``Sig``)  using the default  parameters:
%   embedding dimension = 2, time delay = 1, symbols = 6, logarithm = natural
%   normalisation = false.
%      * Note: ``CondEn(m=1)`` returns the Shannon entropy of ``Sig``.
%
%   [Cond, SEw, SEz] = CondEn(Sig, name, value, ...)
% 
%   Returns the corrected conditional entropy estimates (``Cond``) from the data
%   sequence (``Sig``) using the specified name/value pair arguments:
% 
%      * ``m``     - Embedding Dimension, an integer > 1
%      * ``tau``   - Time Delay, a positive integer
%      * ``c``     - Number of symbols, an integer > 1
%      * ``Logx``  - Logarithm base, a positive scalar 
%      * ``Norm``  - Normalisation of ``Cond`` value, a boolean.
%              * [false] no normalisation - default
%              * [true]  normalises w.r.t Shannon entropy of data sequence ``Sig``  
%
%   See also:
%       XCondEn, MSEn, PermEn, DistEn, XPermEn
%   
%   References:
%     [1] Alberto Porta, et al.,
%           "Measuring regularity by means of a corrected conditional
%           entropy in sympathetic outflow." 
%           Biological cybernetics 
%           78.1 (1998): 71-78.
%           


narginchk(1,11)
Sig = squeeze(Sig);

p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isscalar(x) && (x > 0);
Chk3 = @(x) isscalar(x) && (x > 1) && (mod(x,1)==0);
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'m',2,Chk3);
addParameter(p,'tau',1,Chk);
addParameter(p,'c',6,Chk3);
addParameter(p,'Logx',exp(1),Chk2);
addParameter(p,'Norm',false,@(x) islogical(x));
parse(p,Sig,varargin{:})
m = p.Results.m; tau = p.Results.tau; Logx = p.Results.Logx; 
Norm = p.Results.Norm; c = p.Results.c;

Sig = (Sig-mean(Sig))/std(Sig,1);
Edges = linspace(min(Sig),max(Sig),c+1);
Sx = discretize(Sig,Edges);
N = length(Sx);
SEw = zeros(1,m-1);
SEz = zeros(1,m-1);
Prcm = zeros(1,m-1);
Xi = zeros(N,m);
Xi(:,m) = Sx;
for k = 1:m-1
    Nx = N-(k*tau);
    Xi(1:Nx,end-k) = Sx((k*tau)+1:N);
    Wi = (c.^[k-1:-1:0])*Xi(1:Nx,m-k+1:m)';
    Zi = (c.^[k:-1:0])*Xi(1:Nx,m-k:m)';
    Pw = histcounts(Wi,[min(Wi)-.5:max(Wi)+.5]);
    Pz = histcounts(Zi,[min(Zi)-.5:max(Zi)+.5]);
    Prcm(k) = sum(Pw==1)/Nx;
    
    if sum(Pw)~= Nx || sum(Pz)~= Nx
        warning('Potential error estimating probabilities.')
    end
    
    Pw(Pw==0) = []; Pw = Pw/N;
    Pz(Pz==0) = []; Pz = Pz/N;
    SEw(k) = -Pw*log(Pw)'/log(Logx);
    SEz(k) = -Pz*log(Pz)'/log(Logx);
    clear Pw Pz Wi Zi
end

Temp = histcounts(Sx,c)/N;
Temp(Temp==0) = [];
S1 = -Temp*log(Temp)'/log(Logx);
Cond = SEz - SEw + Prcm*S1;
Cond = [S1 Cond];
if Norm
    Cond = Cond/S1;
end
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub