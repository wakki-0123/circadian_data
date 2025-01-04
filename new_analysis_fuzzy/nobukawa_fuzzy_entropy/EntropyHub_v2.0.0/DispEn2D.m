function [Disp2D, RDE] = DispEn2D(Mat, varargin)
% DispEn2D  estimates the bidimensional dispersion entropy of a data matrix.
%
%   [Disp2D, RDE] = DispEn2D(Mat) 
% 
%   Returns the bidimensional dispersion entropy estimate (``Disp2D``) and 
%   reverse bidimensional dispersion entropy (``RDE``) estimated for the data 
%   matrix (``Mat``) using the default parameters:  time delay = 1, symbols = 3,
%   logarithm = natural, data transform = normalised cumulative density function (ncdf)
%   matrix template size = [floor(H/10) floor(W/10)]  
%   (where H and W represent the height (rows) and width (columns) of the
%   data matrix ``Mat``) 
%   * The minimum number of rows and columns of Mat must be > 10.
%
%   [Disp2D, RDE] = DispEn2D(Mat, name, value, ...)
% 
%   Returns the bidimensional dispersion entropy (``Disp2D``) estimate and
%   reverse bidimensional dispersion entropy (``RDE``) for the data matrix
%   (``Mat``) using the specified name/value pair arguments:
% 
%      * ``m``     - Template submatrix dimensions, an integer scalar 
%        (for submatrix with same height and width) or a two-element vector of
%        integers [height, width] with a minimum value > 1.
%        (default: [floor(H/10) floor(W/10)])
%      * ``tau``   - Time Delay, a positive integer   (default: 1)
%      * ``c``     - Number of symbols, an integer > 1
%      * ``Typex`` - Type of symbolic mapping transform, one of the following:
%        {``linear``, ``kmeans``, ``ncdf``, ``equal``}
%        See the `EntropyHub Guide` for more info on these transforms.
%      * ``Logx``  - Logarithm base, a positive scalar
%      * ``Norm``  - Normalisation of ``Disp2D`` and ``RDE`` values, a boolean:
%                - [false]   no normalisation - default
%                - [true]    normalises w.r.t number of possible dispersion 
%                  patterns (``c^m``).
%      * ``Lock``  - By default, ``DispEn2D`` only permits matrices with a maximum
%        size of 128 x 128 to prevent RAM overload. 
%        e.g. For ``Mat`` = [200 x 200], ``m = 3``, and ``tau = 1``, ``DispEn2D`` 
%        creates a vector of 753049836 elements. To enable matrices
%        greater than [128 x 128] elements, set ``'Lock' = false``. (default: true)
% 
%        **WARNING: unlocking the permitted matrix size may cause memory
%        errors that could lead Matlab to crash.**
%
%   See also:
%       DispEn, SampEn2D, FuzzEn2D, DistEn2D
%   
%   References:
%     [1] Hamed Azami, et al.,
%           "Two-dimensional dispersion entropy: An information-theoretic 
%           method for irregularity analysis of images."
%           Signal Processing: Image Communication, 
%           75 (2019): 178-187.
% 

narginchk(1,15)
Mat = squeeze(Mat);
[NL, NW] = size(Mat);

q = inputParser;
Chk = @(x) isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isvector(x) && (min(x) > 1) && (max(mod(x,1))==0);
Chk3 = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
Chk4 = {'linear','kmeans','ncdf','equal'};
Chkx = @(x) isnumeric(x) && ismatrix(x) && (min(size(x)) > 10);
addRequired(q,'Mat', Chkx);
addParameter(q,'m',[floor(NL/10) floor(NW/10)],Chk2);
addParameter(q,'tau', 1, Chk);
addParameter(q,'c',3,@(x) isnumeric(x) && (x > 1) && (mod(x,1)==0));
addParameter(q,'Logx',exp(1),Chk3);
addParameter(q,'Norm',false,@(x) islogical(x));
addParameter(q,'Typex','ncdf',@(x) ischar(x) && any(validatestring(lower(x),Chk4)));
addParameter(q,'Lock',true,@(x) islogical(x));
parse(q,Mat,varargin{:})
tau = q.Results.tau; Logx = q.Results.Logx;
Norm = q.Results.Norm; Typex = q.Results.Typex;
c = q.Results.c;

if (NL > 128 || NW > 128) && q.Results.Lock
    error(['To prevent memory errors, matrix width & length' ...
        ' must have <= 128 elements. \nTo estimate DistEn2D ', ...
        'for the currect matrix (%dx%d) change Lock to ''unlocked''.', ...
        '\nCaution: unlocking the safe matrix size may cause ', ...
        'MatLab to crash.'],NL,NW)
end
if length(q.Results.m)==1
    mL = q.Results.m; 
    mW = q.Results.m;
else
    mL = q.Results.m(1); 
    mW = q.Results.m(2);
end
if Logx == 0
    Logx = exp(1);
end

switch lower(Typex)
    case 'linear'
        Zi = discretize(Mat,linspace(min(Mat(:)),max(Mat(:)),c+1));
        
    case 'kmeans'         
        
        [Zx,Clux] = kmeans(Mat(:), c, 'MaxIter', 200);
        Zi = zeros(size(Mat(:)));
        [~,xx] = sort(Clux);  
        for k = 1:c
            Zi(Zx==xx(k)) = k;
        end
        
        Zi = reshape(Zi,size(Mat));
        clear Clux Zx xx
        
    case 'ncdf'       
        Zx = normcdf(Mat(:),mean(Mat(:)),std(Mat(:),1));
        Zi = discretize(reshape(Zx,size(Mat)),linspace(0,1,c+1));     
        
%     case 'finesort'
%         Zx = normcdf(Mat(:),mean(Mat(:)),std(Mat(:),1));
%         Zi = discretize(Zx,linspace(0,1,c+1));
%         
%         Ym = zeros(N-(m-1)*tau, m);
%         for n = 1:m
%             Ym(:,n) = Zx(1+(n-1)*tau:N-((m-n)*tau));
%         end
%         Yi = floor(max(abs(diff(Ym,[],2)),[],2)/(rho*std(abs(diff(Sig)),1)));
%         clear Zx Ym
        
    case 'equal'        
        T_Mat = Mat';
        [~,ix] = sort(T_Mat(:));
        xx = round(linspace(0,numel(Mat),c+1),'TieBreaker',"even");        
        Zi = zeros(1,numel(Mat));
        for k = 1:c
            Zi(ix(xx(k)+1:xx(k+1))) = k;
        end
        Zi = reshape(Zi,size(Mat'))';
        clear ix xx
end

NL = NL - (mL - 1)*tau;
NW = NW - (mW - 1)*tau;
X = zeros(NL*NW,mL*mW);
p = 0;
for k = 1:NL
    for n = 1:NW
        p = p+1;
        Temp = Zi(k:tau:(mL-1)*tau+k,n:tau:(mW-1)*tau+n);
        X(p,:) = Temp(:);        
    end
end

if p ~= NL*NW
    warning('Potential error with submatrix division.')
end

T = unique(X,'rows');
Nx = size(T,1);
Counter = zeros(1,Nx);
for n = 1:Nx
    Counter(n) = sum(~any(X - T(n,:),2));
end
Ppi = Counter(Counter~= 0)/size(X,1);
% RDE = sum(Ppi.^2) - (1/Nx);

if c^(mL*mW) > 10^16
    error(sprintf(['RDE cannot be estimated with c = %d and a submatrix of size %d x %d.', ... 
    ' Required floating point precision exceeds 10^16.', ...
    ' Consider reducing the template submatrix size (m) or the number of symbols (c).'],c,mL,mW))
end
    
RDE = sum((Ppi - (1/(c^(mL*mW)))).^2);
if round(sum(Ppi),4) ~= 1
    warning('Potential Error calculating probabilities')
end

Disp2D = -sum(Ppi.*log(Ppi)/log(Logx));
% if Norm
%     Disp2D = Disp2D/(log(Nx)/log(Logx));
%     RDE = RDE/(1 - (1/(Nx)));
% end
if Norm
    Disp2D = Disp2D/(log(c^(mL*mW))/log(Logx));
    RDE = RDE/(1 - (1/(c^(mL*mW))));
end
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub