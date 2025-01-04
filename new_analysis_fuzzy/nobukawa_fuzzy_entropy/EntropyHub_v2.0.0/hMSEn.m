function [MSx, Sn, CI] = hMSEn(Sig, Mobj, varargin) 
% hMSEn  returns the hierarchical entropy of a univariate data sequence.
%
%   [MSx,Sn,CI] = hMSEn(Sig, Mobj) 
% 
%   Returns a vector of entropy values (``MSx``) calculated at each node in the
%   hierarchical tree, the average entropy value across all nodes at each 
%   scale (``Sn``), and the complexity index (``CI``) of the hierarchical tree
%   (i.e. ``sum(Sn)``) for the data sequence (``Sig``) using the parameters specified
%   by the multiscale object (``Mobj``) over 3 temporal scales (default).
%   The entropy values in ``MSx`` are ordered from the root node (S_00) to the
%   Nth subnode at scale T (S_TN): i.e. S_00, S_10, S_11, S_20, S_21, S_22,
%   S_23, S_30, S_31, S_32, S_33, S_34, S_35, S_36, S_37, S_40, ... , S_TN.
%   The average entropy values in ``Sn`` are ordered in the same way, with the
%   value of the root node given first: i.e. S0, S1, S2, ..., ST
%    
%   [MSx,Sn,CI] = hMSEn(Sig, Mobj, name, value, ...)
% 
%   Returns a vector of entropy values (``MSx``) calculated at each node in the
%   hierarchical tree, the average entropy value across all nodes at each 
%   scale (``Sn``), and the complexity index (``CI``) of the entire hierarchical
%   tree for the data sequence (``Sig``) using the following name/value pair 
%   arguments:
%       * ``Scales``   - Number of temporal scales, an integer > 1   (default = 3)
%         At each scale (T), entropy is estimated for 2^(T-1) nodes.
%       * ``RadNew``   - Radius rescaling method, an integer in the range [1 4].
%         When the entropy specified by ``Mobj`` is ``SampEn`` or ``ApEn``, 
%         RadNew rescales the radius threshold in each sub-sequence
%         at each time scale (Xt). If a radius value (``r``) is specified 
%         by ``Mobj``, this becomes the rescaling coefficient, otherwise
%         it is set to 0.2 (default). The value of ``RadNew`` specifies
%         one of the following methods:
%                   * [1]    Standard Deviation          - r*std(Xt)
%                   * [2]    Variance                    - r*var(Xt)
%                   * [3]    Mean Absolute Deviation     - r*mad(Xt)
%                   * [4]    Median Absolute Deviation   - r*mad(Xt,1)
%       * ``Plotx``  - When ``Plotx == true``, returns a plot of the average entropy 
%         value at each time scale (i.e. the multiscale entropy curve)
%         and a hierarchical graph showing the entropy value of each node
%         in the hierarchical tree decomposition.  (default: false)
% 
%   See also:
%       MSobject, MSEn, rMSEn, cMSEn, XMSEn, hXMSEn, rXMSEn, cXMSEn
% 
%   References:
%       [1] Ying Jiang, C-K. Peng and Yuesheng Xu,
%           "Hierarchical entropy analysis for biological signals."
%           Journal of Computational and Applied Mathematics
%           236.5 (2011): 728-742.
% 
narginchk(2,8)
Sig = squeeze(Sig);
p = inputParser;
addRequired(p,'Sig',@(x) isnumeric(x) && isvector(x) && (length(x) > 10));
addRequired(p,'Mobj',@(x) isstruct(x));
addParameter(p,'Scales',3,@(x) isnumeric(x) && (length(x)==1) && (x>1) && mod(x,1)==0);
addParameter(p,'Plotx',false,@(x) islogical(x));
addParameter(p,'RadNew',0,@(x) x==0 || (ismember(x,1:4) && ...
    any(validatestring(func2str(Mobj.Func),{'SampEn';'ApEn'}))));
parse(p,Sig, Mobj, varargin{:})
RadNew = p.Results.RadNew;
Scales = p.Results.Scales;

if strcmp(func2str(Mobj.Func),'SampEn'),  Mobj.Vcp = false; end

Fields = fieldnames(Mobj);
Y = struct2cell(Mobj);
C = [Fields(2:end),Y(2:end)].';

assert(~startsWith(lower(func2str(Mobj.Func)), 'x'),...
    "Base entropy estimator is a cross-entropy method. To perform heirarchical multiscale CROSS-entropy estimation, use hXMSEn.")

if RadNew    
    switch RadNew
        case 1
            Rnew = @(x) std(x,1);
        case 2
            Rnew = @(x) var(x,1);
        case 3            
            Rnew = @(x) mad(x);
        case 4
            Rnew = @(x) mad(x,1);
    end
    
    try
        C_Loc = find(strcmp(C(1,:),'r'));
        Cx = C{2,C_Loc};
    catch
        Cy = {'Standard Deviation';'Variance';...
            'Mean Abs Deviation';'Median Abs Deviation'};
        warning(['No radius value provided.\n' ...
            'Default set to 0.2*(%s) of each new time-series.'],Cy{RadNew})
        C_Loc = size(C,1) + 1;
        %  C_Loc = length(C(1,:)) + 1;
        C{1,C_Loc} = 'r';
        Cx = .2;
    end
end
[XX,N] = Hierarchy(Sig,Scales);
MSx = zeros(1,size(XX,1));
for T = 1:size(XX,1)   
    fprintf(' .')
    Temp = XX(T,1:N/(2^(floor(log2(T)))));
    if RadNew
        C{2,C_Loc} = Cx*Rnew(Temp); 
    end
    Temp2 = Mobj.Func(Temp,C{:});
    MSx(T) = Temp2(end);   
    clear Temp Temp2    
end
Sn = zeros(1,p.Results.Scales);
for t = 1:Scales
    Sn(t) = mean(MSx(2^(t-1):(2^t)-1));
end
CI = sum(Sn);
fprintf(' .\n')
if any(isnan(MSx))
    fprintf('\nSome entropy values may be undefined.')
end

if p.Results.Plotx
    figure, subplot(4,1,1), hold on
    plot(1:Scales,Sn,'color',[8 63 77]/255,'LineWidth',3)
    scatter(1:Scales,Sn,60,[1 0 1],'filled')
    xlabel('Scale Factor','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
    ylabel('Entropy Value','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
    title(sprintf('Hierarchical Multiscale %s Entropy',func2str(Y{1})),...
        'FontSize',16,'FontWeight','bold','Color',[7 54 66]/255)
        
    subplot(4,1,2:4), hold on,
    nodes = [0 repelem(1:(2^(Scales-1))-1,2)];
    [x,y] = treelayout(nodes);
    Tx = 2.^(0:Scales - 1);
    Temp = zeros(Scales, 2^(Scales-1));
    for k = 1:Scales
        Temp(k,:) = Tx(k) + repelem(0:Tx(k)-1,1,Tx(end-k+1));   
    end
    y_new = zeros(size(nodes));
    for i=1:length(nodes)
        [row,~] = find(Temp==i, 1);
        y_new(i) = max(y) - (row-1)*min(y);
    end
    for c = 1:size(Temp, 2)
        line_x = x(Temp(Temp(:,c)>0, c));
        line_y = y_new(Temp(Temp(:,c)>0, c));
        line(line_x, line_y,'Color',[8 63 77]/255,'LineWidth',2.5);
    end
    for i = 1:length(nodes)        
        scatter(x(i),y_new(i),100*(MSx(i)-min(MSx)+1)/(abs(min(MSx))+1),...
            'MarkerEdgeColor',[1 0 1],'MarkerFaceColor',[1 0 1])
        text(x(i),y_new(i),num2str(round(MSx(i),3)))
    end
    
%     MM = zeros(2*p.Results.Scales - 1,T);
%     mid = ceil(T/2);
%     for t = 1:p.Results.Scales
%         MM(2*t -1,mid-(2^(t-1))+1:2:mid-1+(2^(t-1))) = MSx(2^(t-1):(2^t)-1);
%     end  
%     
%     figure, b = bar3(MM);
%     zlim([min(MSx(MSx>0)) max(MSx)]); 
%     cmap = colormap('cool');  colormap([1 1 1; cmap])
%     for k = 1:length(b)
%         zdata = b(k).ZData;
%         b(k).CData = zdata;
%         b(k).FaceColor = 'interp';
%     end
%     axis square, view([30 55])        
%     yticks(1:2:(2*T -1));  yticklabels(0:T-1)
%     xticks(''); %     xticks(1:2:(2*p.Results.Scales - 1))    
%     xlabel('Subtree Nodes','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
%     ylabel('Scale Factor','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
%     zlabel('Entropy','FontSize',12,'FontWeight','bold','Color',[7 54 66]/255)
%     title(sprintf('Hierarchical Multiscale (%s) Entropy',func2str(Y{1})),...
%         'FontSize',16,'FontWeight','bold','Color',[7 54 66]/255)    
end
end

function [U, N] = Hierarchy(Z,sx)
N = 2^floor(log2(length(Z)));
if mod(log2(length(Z)),1)~= 0
    fprintf(['Only first %d samples (2^Scales) were used in hierarchical decomposition.'...
        '\nThe last %d samples of the data sequence were ignored.'], N, length(Z)-N)
end
if N/(2^(sx-1)) < 8
   error(['Data length (%d) is too short to estimate entropy at the lowest '...
       'subtree.\nConsider reducing the number of scales.'],N) 
end

Z = Z(1:N);
U = zeros((2^sx)-1,N);
U(1,:) = Z;
p=2;
for k = 1:sx-1
    for n = 1:2^(k-1)
        Temp = U(2^(k-1)+n-1,:);        
        U(p,1:N/2)  = (Temp(1:2:end) + Temp(2:2:end))/2;
        U(p+1,1:N/2)= (Temp(1:2:end) - Temp(2:2:end))/2;
        p=p+2;
        clear Temp
    end
end
end


%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub
