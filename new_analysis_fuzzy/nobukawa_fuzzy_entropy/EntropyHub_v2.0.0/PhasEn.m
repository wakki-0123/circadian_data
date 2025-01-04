function [Phas] = PhasEn(Sig, varargin)
% PhasEn  estimates the phase entropy of a univariate data sequence.
%
%   [Phas] = PhasEn(Sig) 
% 
%   Returns the phase entropy (``Phas``) estimate of the data sequence (``Sig``)
%   using the default parameters: 
%   angular partitions = 4, time delay = 1, logarithm = natural,
%   normalisation = true
%
%   [Phas] = PhasEn(Sig, name, value, ...)
% 
%   Returns the phase entropy (``Phas``) estimate of the data sequence (``Sig``)  
%   using the specified name/value pair arguments:
%       * ``K``     - Angular partitions (coarse graining), an integer > 1
%                * Note: Division of partitions begins along the positive
%                        x-axis. As this point is somewhat arbitrary, it is
%                        recommended to use even-numbered (preferably
%                        multiples of 4) partitions for sake of symmetry. 
%       * ``tau``   - Time Delay, a positive integer
%       * ``Logx``  - Logarithm base, a positive scalar  
%       * ``Norm``  - Normalisation of Phas value, a boolean:
%                * [false] no normalisation
%                * [true]  normalises w.r.t. the # partitions (``Log(K)``) (Default)
%       * ``Plotx`` - When ``Plotx == true``, returns Poincaré plot (default: false)
%
%   See also:
%       SampEn, ApEn, GridEn, MSEn, SlopEn, CoSiEn, BubbEn
%   
%   References:
%     [1] Ashish Rohila and Ambalika Sharma,
%           "Phase entropy: a new complexity measure for heart rate
%           variability." 
%           Physiological measurement
%           40.10 (2019): 105006.
% 


narginchk(1,11)
Sig = squeeze(Sig);

p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
Chk3 = @(x) isnumeric(x) && isscalar(x) && (x > 1) && (mod(x,1)==0);
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'K',4,Chk3);
addParameter(p,'tau',1,Chk);
addParameter(p,'Logx',exp(1),Chk2);
addParameter(p,'Norm',true,@(x) islogical(x));
addParameter(p,'Plotx',false,@(x) islogical(x));
parse(p,Sig,varargin{:})
K = p.Results.K; tau = p.Results.tau; Logx = p.Results.Logx; 
Norm = p.Results.Norm; Plotx = p.Results.Plotx; 

if size(Sig,1) < size(Sig,2)
    Sig = Sig';
end
Yn = Sig(1+2*tau:end) - Sig(tau+1:end-tau);
Xn = Sig(tau+1:end-tau) - Sig(1:end-2*tau);
Theta_r = atan(Yn./Xn);
Theta_r(Yn<0 & Xn<0) = Theta_r(Yn<0 & Xn<0) + pi;
Theta_r(Yn<0 & Xn>0) = Theta_r(Yn<0 & Xn>0) + 2*pi;
Theta_r(Yn>0 & Xn<0) = Theta_r(Yn>0 & Xn<0) + pi;

Limx = ceil(max(abs([Yn; Xn])));
Angs = linspace(0,2*pi,K+1);
Tx = zeros(K,length(Theta_r));
Si = zeros(1,K);
for n = 1:K
    Temp = (Theta_r > Angs(n) & Theta_r < Angs(n+1));
    Tx(n,Temp) = 1;
    Si(n) = sum(Theta_r(Temp));
end

Si(Si==0) = [];
Phas = -sum((Si/sum(Si)).*(log(Si/sum(Si))/log(Logx)));
if Norm
    Phas = Phas/(log(K)/log(Logx));
end
if Plotx   
    Tx = logical(Tx);
    Ys = sin(Angs)*Limx*sqrt(2);
    Xs = cos(Angs)*Limx*sqrt(2);
    Cols = [zeros(1,K);repmat(randperm(K)/K,2,1)]';    
    figure(), hold on
    for n = 1:K
        plot(Xn(Tx(n,:)),Yn(Tx(n,:)),'.','Color',Cols(n,:))
    end
    plot([zeros(1,K+1);Xs],[zeros(1,K+1);Ys],'m')
    axis([-Limx Limx -Limx Limx],'square')
    xlabel('X(n + \tau)  -  X(n)'), ylabel('X(n + 2\tau)  -  X(n + \tau)')
    xticks([-Limx 0 Limx]), yticks([-Limx 0 Limx])    
end
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub