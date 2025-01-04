function [MSx, CI] = MSEn2D(Mat, Mobj, varargin) 
% MSEn2D  returns the multiscale entropy of a bidimensional data matrix.
%
%   [MSx,CI] = MSEn2D(Mat, Mobj) 
% 
%   Returns a vector of multiscale entropy values (``MSx``) and the complexity 
%   index (``CI``) of the data matrix ``Mat`` using the parameters specified 
%   by the multiscale object (``Mobj``) over 3 temporal scales with coarse-
%   graining (default). 
% 
%   [MSx,CI] = MSEn2D(Mat, Mobj, name, value, ...)
% 
%   Returns a vector of multiscale entropy values (``MSx``) and the complexity 
%   index (``CI``) of the data matrix ``Mat`` using the parameters specified by
%   the multiscale object (``Mobj``) and the following name/value pair arguments:
%   
%       * ``Scales``   - Number of temporal scales, an integer > 1   (default = 3)
%       * ``Methodx``  - Graining method, one of the following: [default = ``'coarse'``]
%         {``'coarse'``, ``'modified'``} 
%       * ``RadNew``   - Radius threshold rescaling method, an integer in the range [1 4].
%         When the entropy specified by ``Mobj`` is ``SampEn2D``, 
%         RadNew rescales the radius threshold in each sub-matrix
%         at each time scale (Xt). If a radius value (``r``) is specified 
%         by ``Mobj``, this becomes the rescaling coefficient, otherwise
%         it is set to 0.2 (default). The value of ``RadNew`` specifies
%         one of the following methods:
%                   * [1]    Standard Deviation          - r*std(Xt)
%                   * [2]    Variance                    - r*var(Xt)
%                   * [3]    Mean Absolute Deviation     - r*mad(Xt)
%                   * [4]    Median Absolute Deviation   - r*mad(Xt,1)
%       * ``Plotx``    - When ``Plotx == true``, returns a plot of the entropy value at
%         each time scale (i.e. the multiscale entropy curve) [default: false]
% 
%   For further information on multiscale graining procedures, see
%   the `EntropyHub guide <https://github.com/MattWillFlood/EntropyHub/blob/main/EntropyHub%20Guide.pdf>`_
% 
%   See also:
%       MSobject, SampEn2D, DispEn2D, DistEn2D, FuzzEn2D, EspEn2D, PermEn2D
% 
%   References:
%      [1] Luiz E.V. Silva, et al.,
%           "Two-dimensional multiscale entropy analysis: Applications to 
%           image texture evaluation."
%           Signal Processing
%           147 (2018): 224-232.
% 
%      [2] Cristina Morel and Anne Humeau-Heurtier
%           "Multiscale permutation entropy for two-dimensional patterns." 
%           Pattern Recognition Letters 
%           150 (2021): 139-146.
% 
%      [3] Mirvana Hilal and Anne Humeau-Heurtier
%           "Bidimensional Fuzzy Entropy: Principle Analysis and 
%           Biomedical Applications"
%           IEEE Engineering in Medicine and Biology Society (EMBS) Conference
%           2019: 4811-4814
% 

narginchk(2,10)
Mat = squeeze(Mat);
p = inputParser;
Chk = {'coarse';'modified'};
addRequired(p,'Mat',@(x) isnumeric(x) && (min(size(x))>10));
addRequired(p,'Mobj',@(x) isstruct(x));
addParameter(p,'Methodx','coarse',@(x) any(validatestring(string(lower(x)),Chk)));
addParameter(p,'Scales',3,@(x) isnumeric(x) && (length(x)==1) && (x>1));
addParameter(p,'Plotx',false,@(x) islogical(x));
addParameter(p,'RadNew',0,@(x) ismember(x,0:4) && ...
    any(validatestring(func2str(Mobj.Func),{'SampEn2D'})));
parse(p,Mat, Mobj, varargin{:})
Methodx = str2func(lower(p.Results.Methodx));
MSx = zeros(1,p.Results.Scales);
RadNew = p.Results.RadNew;
Fields = fieldnames(Mobj);
Y = struct2cell(Mobj);
C = [Fields(2:end),Y(2:end)].';

if RadNew    
    switch RadNew
        case 1
            Rnew = @(x) std(x(:),1);
        case 2
            Rnew = @(x) var(x(:),1);
        case 3            
            Rnew = @(x) mad(x(:));
        case 4
            Rnew = @(x) mad(x(:),1);
    end
    
    try
        C_Loc = find(strcmp(C(1,:),'r'));
        Cx = C{2,C_Loc};
    catch
        Cy = {'Standard Deviation';'Variance';...
            'Mean Abs Deviation';'Median Abs Deviation'};
        warning(['No radius value provided.\n' ...
            'Default set to 0.2*(%s) of each new time-series.'],Cy{RadNew})
        %  C_Loc = length(C(1,:)) + 1;
        C_Loc = size(C,1) + 1;
        
        C{1,C_Loc} = 'r';
        Cx = .2;
    end
end
for T = 1:p.Results.Scales       
    fprintf(' .')
    Temp = Methodx(Mat,T);  
    if RadNew
        C{2,C_Loc} = Cx*Rnew(Temp);
    end
    Temp2 = Mobj.Func(Temp,C{:});
    MSx(T) = Temp2(end);  
    clear Temp Temp2    
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
       func2str(Methodx)),'FontSize',16,'FontWeight','bold','Color',[7 54 66]/255)
   ylim([0 max(MSx)+.2])
end
end

function [Y] = coarse(Z,sx)
    Nh = floor(size(Z,1)/sx);
    Nw = floor(size(Z,2)/sx);
    Y = squeeze(mean(mean(reshape(Z(1:sx*Nh,1:sx*Nw),sx,Nh,sx,Nw),1),3));
end
function [Y] = modified(Z,sx)   
    Y = conv2(Z, ones(sx), 'valid')/(sx*sx);
end


%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub