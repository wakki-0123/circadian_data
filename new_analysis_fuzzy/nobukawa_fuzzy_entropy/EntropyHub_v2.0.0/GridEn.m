function [GDE, GDR, PIx, GIx, SIx, AIx] = GridEn(Sig, varargin) 
% GridEn  estimates the gridded distribution entropy of a univariate data sequence.
%
%   [GDE, GDR] = GridEn(Sig) 
% 
%   Returns the gridded distribution entropy (``GDE``) and the gridded 
%   distribution rate (``GDR``) estimated from the data sequence (``Sig``) using 
%   the default  parameters:
%   grid coarse-grain = 3, time delay = 1, logarithm = natural
%  
%   [GDE, GDR, PIx, GIx, SIx, AIx] = GridEn(Sig)
% 
%   In addition to ``GDE`` and ``GDR``, ``GridEn`` returns the following indices 
%   estimated from the data sequence (``Sig``) using the default  parameters:
%     -  ``PIx``   - Percentage of points below the line of identity (LI)
%     -  ``GIx``   - Proportion of point distances above the LI
%     -  ``SIx``   - Ratio of phase angles (w.r.t. LI) of the points above the LI
%     -  ``AIx``   - Ratio of the cumulative area of sectors of points above the LI
%
%   [GDE, GDR, ...,] = GridEn(Sig, name, value, ...)
% 
%   Returns the gridded distribution entropy (``GDE``) estimate of the data 
%   sequence (``Sig``) using the specified name/value pair arguments:
% 
%      * ``m``     - Grid coarse-grain (``m`` x ``m`` sectors), an integer > 1
%      * ``tau``   - Time Delay, a positive integer
%      * ``Logx``  - Logarithm base, a positive scalar
%      * ``Plotx`` - When ``Plotx == true``, returns gridded Poicare plot and
%        a bivariate histogram of the grid point distribution (default: false)  
%
%   See also:
%       PhasEn, CoSiEn, SlopEn, BubbEn, MSEn
%   
%   References:
%     [1] Chang Yan, et al.,
%           "Novel gridded descriptors of Poincare plot for analyzing 
%           heartbeat interval time-series." 
%           Computers in biology and medicine 
%           109 (2019): 280-289.
% 
%     [2] Chang Yan, et al. 
%           "Area asymmetry of heart rate variability signal." 
%           Biomedical engineering online 
%           16.1 (2017): 1-14.
% 
%     [3] Alberto Porta, et al.,
%           "Temporal asymmetries of short-term heart period variability 
%           are linked to autonomic regulation." 
%           American Journal of Physiology-Regulatory, Integrative and 
%           Comparative Physiology 
%           295.2 (2008): R550-R557.
% 
%     [4] C.K. Karmakar, A.H. Khandoker and M. Palaniswami,
%           "Phase asymmetry of heart rate variability signal." 
%           Physiological measurement 
%           36.2 (2015): 303.
% 


narginchk(1,9)
Sig = squeeze(Sig);
if size(Sig,1) > size(Sig,2)
    Sig = Sig';
end

p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
Chk3 = @(x) isnumeric(x) && isscalar(x) && (x > 1) && (mod(x,1)==0);
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addParameter(p,'m',3,Chk3);
addParameter(p,'tau',1,Chk);
addParameter(p,'Logx',exp(1),Chk2);
addParameter(p,'Plotx',false,@(x) islogical(x));
parse(p,Sig,varargin{:})
m = p.Results.m; tau = p.Results.tau; 
Logx = p.Results.Logx; Plotx = p.Results.Plotx; 

if Logx == 0
    Logx = exp(1);
end
Sig_n = (Sig-min(Sig))/range(Sig);
Temp = [Sig_n(1:end-tau);Sig_n(1+tau:end)];
N = hist3(Temp',[m,m]);
Pj = flipud(N')/size(Temp,2); Ppi = Pj(Pj>0);
if round(sum(Ppi)) ~= 1
    warning('Potential error of estimated probabilities: P = %d', sum(Ppi))
end
GDE = -sum(Ppi.*(log(Ppi)/log(Logx)));
GDR = sum(N(:)~=0)/(m*m);

if nargout > 2
    T2   = atand(Temp(2,:)./Temp(1,:))';
    Dup  = sum(abs(diff(Temp(:,T2>45),[],1)));
    Dtot = sum(abs(diff(Temp(:,T2~=45),[],1)));
    Sup  = sum((T2(T2>45)-45));
    Stot = sum(abs(T2(T2~=45)-45));
    Aup  = sum(abs(((T2(T2>45)-45)).*sqrt(sum(Temp(:,T2>45).^2)')));
    Atot = sum(abs(((T2(T2~=45)-45)).*sqrt(sum(Temp(:,T2~=45).^2)'))); 
    PIx = 100*sum(T2 < 45)/sum(T2~=45);
    GIx = 100*Dup/Dtot;
    SIx = 100*Sup/Stot;
    AIx = 100*Aup/Atot;
end
if Plotx
    ntrvl = linspace(0,1,m+1);
    figure, subplot(1,2,1), hold on,
    plot(Sig_n(1:end-tau),Sig_n(tau+1:end),'.m')
    plot([ntrvl;ntrvl],[zeros(1,m+1);ones(1,m+1)],'c')
    plot([zeros(1,m+1);ones(1,m+1)],[ntrvl;ntrvl],'c')
    plot([0 1],[0 1],'k'),    axis square
    xlabel('X_n'), ylabel('X_n _+ _\tau'), xticks([0 1]), yticks([0 1])

    subplot(1,2,2)
    histogram2(Temp(1,:),Temp(2,:),m, 'DisplayStyle','tile','ShowEmptyBins','on')
    colormap('cool'), axis square
    xlabel('X_n'), ylabel('X_n _+ _\tau'), xticks([0 1]), yticks([0 1])
end
end
%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub
