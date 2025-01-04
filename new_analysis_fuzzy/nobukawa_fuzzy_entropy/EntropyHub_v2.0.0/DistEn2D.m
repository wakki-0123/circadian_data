function [Dist2D] = DistEn2D(Mat, varargin)
% DistEn2D  estimates the bidimensional distribution entropy of a data matrix.
%
%   [Dist2D] = DistEn2D(Mat) 
% 
%   Returns the bidimensional distribution entropy estimate (``Dist2D``)
%   estimated for the data matrix (``Mat``) using the default parameters:
%   time delay = 1, binning method = ``'sturges'``, logarithm = natural, 
%   matrix template size = [floor(H/10) floor(W/10)]  
%   (where H and W represent the height (rows) and width (columns) of the
%   data matrix ``Mat``) 
%   * The minimum number of rows and columns of ``Mat`` must be > 10.
%
%   [Dist2D] = DistEn2D(Mat, name, value, ...)
% 
%   Returns the bidimensional distribution entropy (``Dist2D``) estimate for 
%   the data matrix (``Mat``) using the specified name/value pair arguments:
% 
%      * ``m``     - Template submatrix dimensions, an integer scalar (for sub-
%        matrix with same height and width) or a two-element vector of
%        integers [height, width] with a minimum value > 1.
%        (default: [floor(H/10) floor(W/10)])
%      * ``tau``   - Time Delay, a positive integer   (default: 1)
%      * ``Bins``  - Histogram bin selection method for distance distribution,
%        an integer > 1 indicating the number of bins, or one of the following strings
%        {``'sturges'``, ``'sqrt'``, ``'rice'``, ``'doanes'``} [default: ``'sturges'``]
%      * ``Logx``  - Logarithm base, a positive scalar    (default: natural)
%      * ``Norm``  - Normalisation of ``Dist2D`` value, one of the following integers:
%               -  [0]  no normalisation.
%               -  [1]  normalises values of data matrix (``Mat``) to range [0 1].
%               -  [2]  normalises values of data matrix (``Mat``) to range [0 1],
%                  and normalises the distribution entropy value (``Dist2D``)
%                  w.r.t the number of histogram bins.  [default]
%               -  [3]  normalises the bidimensional distribution entropy value 
%                  (``Dist2D``) w.r.t the number of histogram bins, without 
%                  normalising data matrix values.
%      * ``Lock``  - By default, ``DistEn2D`` only permits matrices with a maximum
%        size of 128 x 128 to prevent RAM overload. 
%        e.g. For ``Mat`` = [200 x 200], ``m = 3``, and ``tau = 1``, ``DistEn2D`` 
%        creates a vector of 753049836 elements. To enable matrices
%        greater than [128 x 128] elements, set ``'Lock' = false``. (default: true)
% 
%        **WARNING: unlocking the permitted matrix size may cause memory
%        errors that could lead Matlab to crash.**
%
%   See also:
%       DistEn, XDistEn, SampEn2D, FuzzEn2D, MSEn
%   
%   References:
%     [1] Hamed Azami, Javier Escudero and Anne Humeau-Heurtier,
%           "Bidimensional distribution entropy to analyze the irregularity
%           of small-sized textures."
%           IEEE Signal Processing Letters 
%           24.9 (2017): 1338-1342.
% 


narginchk(1,13)
Mat = squeeze(Mat);
[NL, NW] = size(Mat);

q = inputParser;
Chk = @(x) isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isvector(x) && (min(x) > 1) && (max(mod(x,1))==0);
Chk3 = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
Chk4 = {'sturges','sqrt','rice','doanes'};
Chk5 = @(x) (isnumeric(x) && (length(x)==1) && (mod(x,1)==0)) ...
    || (ischar(x) && any(validatestring(lower(x),Chk4)));
Chkx = @(x) isnumeric(x) && ismatrix(x) && (min(size(x)) > 10);
addRequired(q,'Mat', Chkx);
addParameter(q,'m',[floor(NL/10) floor(NW/10)],Chk2);
addParameter(q,'tau',1, Chk);
addParameter(q,'Bins','sturges',Chk5);
addParameter(q,'Norm',2,@(x) ismember(x,[0:3]));
addParameter(q,'Logx',2,Chk3);
addParameter(q,'Lock',true,@(x) islogical(x));
parse(q,Mat,varargin{:})
tau = q.Results.tau; Logx = q.Results.Logx;
Bins = q.Results.Bins; Norm = q.Results.Norm;

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
if Norm == 1 || Norm == 2
    Mat = (Mat - min(Mat(:)))/range(Mat(:));
end
if Logx == 0
    Logx = exp(1);
end

NL = NL - (mL - 1)*tau;
NW = NW - (mW - 1)*tau;
X = zeros(NL*NW,mL,mW);
p = 0;
for k = 1:NL
    for n = 1:NW
        p = p+1;
        X(p,:,:) = Mat(k:tau:(mL-1)*tau+k,n:tau:(mW-1)*tau+n);
    end
end

if p ~= NL*NW
    warning('Potential error with submatrix division.')
end
Ny = p*(p-1)/2;
if Ny > 300000000
    warning('Number of pairwise distance calculations is %d', Ny)
end

Y = zeros(1,Ny);
for k = 1:p-1
    Ix = [((k-1)*(p - k/2)+1), k*(p-((k+1)/2))];
    Y(Ix(1):Ix(2)) = max(abs(X(k+1:end,:,:) - X(k,:,:)),[],[2,3]);
end

if ischar(Bins)
    switch lower(Bins)
        case 'sturges'
            Bx = ceil(log2(Ny) + 1);
        case 'rice'
            Bx = ceil(2*(Ny^(1/3)));
        case 'sqrt'
            Bx = ceil(sqrt(Ny));
        case 'doanes'
            sigma = sqrt(6*(Ny-2)/((Ny+1)*(Ny+3)));
            Bx = ceil(1+log2(Ny)+log2(1+abs(skewness(Y)/sigma)));
        otherwise
            error('Please enter a valid binning method')
    end
else
    Bx = Bins;
end
By = linspace(min(Y),max(Y),Bx+1);
Ppi = histcounts(Y,By)/Ny;

if round(sum(Ppi),6) ~= 1
    warning('Potential error estimating probabilities (p=%d).', sum(Ppi))
    Ppi(Ppi==0)=[];
elseif any(Ppi==0)
    fprintf('Note: %d/%d bins were empty',sum(Ppi==0),numel(Ppi));
    Ppi(Ppi==0)=[];
end
Dist2D = -sum(Ppi.*log(Ppi)/log(Logx));
if Norm >= 2
    Dist2D = Dist2D/(log(Bx)/log(Logx));
end
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub