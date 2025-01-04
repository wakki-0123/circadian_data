function [MSx, CI] = XMSEn(Sig1, Sig2, Mobj, varargin) 
% XMSEn  returns the multiscale cross-entropy between two univariate data sequences.
%
%   [MSx,CI] = XMSEn(Sig1, Sig2, Mobj) 
% 
%   Returns a vector of multiscale cross-entropy values (``MSx``) and the
%   complexity index (``CI``) between the data sequences contained in ``Sig1``
%   and ``Sig2`` using the parameters specified by the multiscale object 
%   (``Mobj``) over 3 temporal scales with coarse-graining (default).
%    
%   [MSx,CI] = XMSEn(Sig1, Sig2, Mobj, name, value, ...)
% 
%   Returns a vector of multiscale cross-entropy values (``MSx``) and the
%   complexity index (``CI``) between the data sequences contained in ``Sig1``
%   and ``Sig2`` using the parameters specified by the multiscale object 
%   (``Mobj``) and the following name/value pair arguments:
% 
%       * ``Scales``   - Number of temporal scales, an integer > 1   (default = 3)
%       * ``Methodx``  - Graining method, one of the following: [default = ``'coarse'``]
%         {``'coarse'``,``'generalized'``, ``'modified'``, ``'imf'``, ``'timeshift'``} 
%       * ``RadNew``   - Radius rescaling method, an integer in the range [1 4].
%         When the entropy specified by ``Mobj`` is ``XSampEn`` or ``XApEn``, 
%         ``RadNew`` rescales the radius threshold in each sub-sequence
%         at each time scale (Xt). If a radius value (``r``) is specified 
%         by ``Mobj``, this becomes the rescaling coefficient, otherwise
%         it is set to 0.2 (default). The value of ``RadNew`` specifies
%         one of the following methods:
%                   * [1]    Pooled Standard Deviation       - r*std(Xt)
%                   * [2]    Pooled Variance                 - r*var(Xt)
%                   * [3]    Total Mean Absolute Deviation   - r*mad(Xt)
%                   * [4]    Total Median Absolute Deviation - r*mad(Xt,1)
%       * ``Plotx``    - When ``Plotx == true``, returns a plot of the entropy value at
%         each time scale (i.e. the multiscale entropy curve) [default: false]
% 
%   See also MSobject, XSampEn, XApEn, rXMSEn, cXMSEn, hXMSEn, MSEn
%   
%   References:
%     [1] Rui Yan, Zhuo Yang, and Tao Zhang,
%           "Multiscale cross entropy: a novel algorithm for analyzing two
%           time series." 
%           5th International Conference on Natural Computation. 
%           Vol. 1, pp: 411-413 IEEE, 2009.
% 
%     [2] Madalena Costa, Ary Goldberger, and C-K. Peng,
%           "Multiscale entropy analysis of complex physiologic time series."
%           Physical review letters
%           89.6 (2002): 068102.
% 
%     [3] Vadim V. Nikulin, and Tom Brismar,
%           "Comment on “Multiscale entropy analysis of complex physiologic
%           time series”." 
%           Physical review letters 
%           92.8 (2004): 089803.
% 
%     [4] Madalena Costa, Ary L. Goldberger, and C-K. Peng. 
%           "Costa, Goldberger, and Peng reply." 
%           Physical Review Letters
%           92.8 (2004): 089804.
% 
%     [5] Antoine Jamin, et al,
%           "A novel multiscale cross-entropy method applied to navigation 
%           data acquired with a bike simulator." 
%           41st annual international conference of the IEEE EMBC
%           IEEE, 2019.
% 
%     [6] Antoine Jamin and Anne Humeau-Heurtier. 
%           "(Multiscale) Cross-Entropy Methods: A Review." 
%           Entropy 
%           22.1 (2020): 45.
%
%     [7] Magdalena Costa and Ary Goldberger, 
%           "Generalized multiscale entropy analysis: Application to
%           quantifying the complex volatility of human heartbeat time
%           series"
%           Entropy
%           17 (2015): 1197–1203
% 

narginchk(3,11)
p = inputParser;
Chk = {'coarse';'modified';'imf';'timeshift';'generalized'};
Chk2 = @(x) isnumeric(x) && isvector(x) && numel(x)>=10;
addRequired(p,'Sig1',Chk2);
addRequired(p,'Sig2',Chk2);
addRequired(p,'Mobj',@(x) isstruct(x));
addParameter(p,'Methodx','coarse',@(x) any(validatestring(string(x),Chk)));
addParameter(p,'Scales',3,@(x) isnumeric(x) && (x>0) && isscalar(x));
addParameter(p,'RadNew',0,@(x) x==0 || (ismember(x,1:4) && ...
    any(validatestring(func2str(Mobj.Func),{'XSampEn';'XApEn'}))));
addParameter(p,'Plotx',false,@(x) islogical(x));
parse(p,Sig1, Sig2, Mobj, varargin{:})
Methodx = str2func(p.Results.Methodx);
RadNew = p.Results.RadNew;
MSx = zeros(1,p.Results.Scales);

if strcmp(func2str(Mobj.Func),'XSampEn'),  Mobj.Vcp = false; end

Fields = fieldnames(Mobj);
Y = struct2cell(Mobj);
C = [Fields(2:end),Y(2:end)].';
Sig1 = Sig1(:);  Sig2 = Sig2(:);

if strcmp(p.Results.Methodx,'imf')
    [Imfx,Resx] = emd(Sig1,'MaxNumIMF',p.Results.Scales-1,'Display',0); 
    [Imfy,Resy] = emd(Sig2,'MaxNumIMF',p.Results.Scales-1,'Display',0); 
    % Sig = zeros(size(Sig,1),2,p.Results.Scales);
    Sig1 = [Imfx Resx];
    Sig2 = [Imfy Resy];    
    clear Imfx Resx Imfy Resy
end

if RadNew    
    switch RadNew
        case 1
            Rnew = @(x,y) sqrt((var(x,1)*(length(x)-1) + var(y,1)*(length(y)-1))/(length(x)+length(y)-1));
        case 2
            Rnew = @(x,y) (var(x,1)*(length(x)-1) + var(y,1)*(length(y)-1))/(length(x)+length(y)-1);
        case 3            
            Rnew = @(x,y) mad([x(:); y(:)]);
        case 4
            Rnew = @(x,y) mad([x(:); y(:)],1);
    end
    
    try
        C_Loc = find(strcmp(C(1,:),'r'));
        Cx = C{2,C_Loc};
    catch
        Cy = {'Standard Deviation';'Variance';...
            'Mean Abs Deviation';'Median Abs Deviation'};
        warning(['No radius value provided.\n' ...
            'Default set to 0.2*(%s) of each new time-series.'],Cy{RadNew})
        %C_Loc = length(C(1,:)) + 1;
        C_Loc = size(C,1) + 1;
        C{1,C_Loc} = 'r';
        Cx = .2;
    end
end

for T = 1:p.Results.Scales    
    fprintf(' .')
    [TempA, TempB] = Methodx(Sig1, Sig2, T);      
    if strcmpi(p.Results.Methodx,'timeshift')
        Tempx = zeros(1,T);
        for k = 1:T
            fprintf(' .')
            if RadNew
                C{2,C_Loc} = Cx*Rnew(TempA(k,:),TempB(k,:));
            end
            Tempy = Mobj.Func(TempA(k,:), TempB(k,:), C{:});  
            Tempx(k) = Tempy(end);
        end
        Temp2 = mean(Tempx);
        clear Tempx Tempy
    else
        if RadNew
            C{2,C_Loc} = Cx*Rnew(TempA, TempB);
        end
        Temp2 = Mobj.Func(TempA, TempB, C{:});
    end
    MSx(T) = Temp2(end);   
    clear TempA TempB Temp2
end
CI = sum(MSx);
fprintf(' .\n')
if any(isnan(MSx))
    fprintf('\nSome entropy values may be undefined.')
end

if p.Results.Plotx
   figure, hold on
   plot(1:p.Results.Scales,MSx,'color',[8 63 77]/255,'LineWidth',3)   
   scatter(1:p.Results.Scales,MSx,60,[1 0 1],'filled')   
   xlabel('Scale Factor','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
   ylabel('Entropy Value','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
   title(sprintf('Multiscale %s (%s-graining method)',func2str(Y{1}),...
       p.Results.Methodx),'FontSize',16,'FontWeight','bold','Color',[7 54 66]/255)
end
end

function [Y1,Y2] = coarse(Za, Zb, sx)
    Ns = floor(length(Za)/sx);
    Y1 = mean(reshape(Za(1:sx*Ns),sx,Ns),1);
    Ns = floor(length(Zb)/sx);
    Y2 = mean(reshape(Zb(1:sx*Ns),sx,Ns),1);
end
function [Y1, Y2] = modified(Za, Zb, sx)
    % Y = movmean(Z,sx);
    % Y = Y(ceil((sx+1)/2):end-ceil(sx/2)+1,:);
    Y1 = conv(Za,ones(1,sx),'valid')/sx;
    Y2 = conv(Zb,ones(1,sx),'valid')/sx;
end
function [Y1, Y2] = imf(Za, Zb, sx)
    Y1 = squeeze(sum(Za(:,1:sx),2));    
    Y2 = squeeze(sum(Zb(:,1:sx),2));
end
function [Y1, Y2] = timeshift(Za, Zb, sx)
    % Y1 = zeros(sx,floor(length(Za)/sx),2);
    Y1 = reshape(Za(1:sx*floor(length(Za)/sx)),sx,floor(length(Za)/sx));
    Y2 = reshape(Zb(1:sx*floor(length(Zb)/sx)),sx,floor(length(Zb)/sx));    
end
function [Y1, Y2] = generalized(Za, Zb, sx)
Na = floor(length(Za)/sx);
Nb = floor(length(Zb)/sx);
Y1 = var(reshape(Za(1:sx*Na),sx,Na),1,1);
Y2 = var(reshape(Zb(1:sx*Nb),sx,Nb),1,1);
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub