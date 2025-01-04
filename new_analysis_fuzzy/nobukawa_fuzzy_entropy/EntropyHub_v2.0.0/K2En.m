function [K2, Ci] = K2En(Sig, varargin) 
% K2En  estimates the Kolmogorov (K2) entropy of a univariate data sequence.
%
%   [K2, Ci] = K2En(Sig) 
% 
%   Returns the Kolmogorov entropy estimates (``K2``) and the correlation
%   integrals (``Ci``) for ``m`` = [1,2] estimated from the data sequence (``Sig``)
%   using the default parameters: embedding dimension = 2, time delay = 1, 
%   distance threshold (``r``) = 0.2*SD(``Sig``), logarithm = natural
%
%   [K2, Ci] = K2En(Sig, name, value, ...)
% 
%   Returns the Kolmogorov entropy estimates (``K2``) for dimensions = [1, ..., ``m``]
%   estimated from the data sequence (``Sig``) using the specified name/value pair
%   arguments:
% 
%       * ``m``     - Embedding Dimension, a positive integer
%       * ``tau``   - Time Delay, a positive integer
%       * ``r``     - Radius Distance Threshold, a positive scalar  
%       * ``Logx``  - Logarithm base, a positive scalar  
%
%   See also:
%       DistEn, XK2En, MSEn
%   
%   References:
%     [1] Peter Grassberger and Itamar Procaccia,
%           "Estimation of the Kolmogorov entropy from a chaotic signal." 
%           Physical review A 28.4 (1983): 2591.
% 
%     [2] Lin Gao, Jue Wang  and Longwei Chen
%           "Event-related desynchronization and synchronization 
%           quantification in motor-related EEG by Kolmogorov entropy"
%           J Neural Eng. 2013 Jun;10(3):03602
% 

narginchk(1,9)
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

N   = length(Sig);
m   = m+1;
Zm  = zeros(N,m);
Ci = zeros(1,m);
for n = 1:m
    N2 = N-(n-1)*tau;
    Zm(1:N2,n) = Sig((n-1)*tau + 1:N);   
    
    % Norm = zeros(N2-1);   
    Norm = inf*ones(N2-1); 
    for k = 1:N2-1
        Temp = repmat(Zm(k,1:n),N2-k,1) - Zm(k+1:N2,1:n);
        Norm(k,k:N2-1) = sqrt(sum(Temp.*Temp,2)); 
    end
    % Norm(Norm==0) = inf;
    Ci(n) = 2*sum(sum(Norm < r))/(N2*(N2-1));     
end
 
K2 = (log(Ci(1:m-1)./Ci(2:m))/log(Logx))/tau;
K2(isinf(K2)) = NaN;
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub