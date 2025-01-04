function [Esp2D] = EspEn2D(Mat, varargin) 
% EspEn2D  estimates the bidimensional Espinosa entropy of a data matrix.
%
%   [Esp2D] = EspEn2D(Mat) 
% 
%   Returns the bidimensional Espinosa entropy estimate (``Esp2D``) 
%   estimated for the data matrix (``Mat``) using the default parameters: 
%   time delay = 1, tolerance threshold = 20, percentage similarity = 0.7
%   logarithm = natural, matrix template size = [floor(H/10) floor(W/10)]  
%   (where H and W represent the height (rows) and width (columns) of 
%   the data matrix ``Mat``) 
%   ** The minimum number of rows and columns of ``Mat`` must be > 10.
%
%   [Esp2D] = EspEn2D(Mat, name, value, ...)
% 
%   Returns the bidimensional Espinosa entropy (``Esp2D``) estimates for the data
%   matrix (``Mat``) using the specified name/value pair arguments:
% 
%      * ``m``     - Template submatrix dimensions, an integer scalar (for sub-
%        matrix with same height and width) or a two-element vector of
%        integers [height, width] with a minimum value > 1.
%        (default: [floor(H/10) floor(W/10)])
%      * ``tau``   - Time Delay, a positive integer   (default: 1)
%      * ``r``     - Tolerance Threshold, a positive scalar (default: 20)
%      * ``ps``    - Percentage similarity, a scalar in range [0 1] (default: 0.7) 
%      * ``Logx``  - Logarithm base, a positive scalar    (default: natural)
%      * ``Lock``  - By default, ``SampEn2D`` only permits matrices with a maximum
%        size of 128 x 128 to prevent RAM overload. 
%        e.g. For ``Mat`` = [200 x 200], ``m = 3``, and ``tau = 1``, ``EspEn2D`` 
%        creates a vector of 753049836 elements. To enable matrices
%        greater than [128 x 128] elements, set ``'Lock' = false``. (default: true)
% 
%        **WARNING: unlocking the permitted matrix size may cause memory
%        errors that could lead Matlab to crash.**
%
%   See also:
%       SampEn2D, FuzzEn2D, DistEn2D, DispEn2D.
% 
%   References:
%      [1] Ricardo Espinosa, et al.,
%           "Two-dimensional EspEn: A New Approach to Analyze Image Texture 
%           by Irregularity." 
%           Entropy,
%           23:1261 (2021)
% 

narginchk(1,13)
Mat = squeeze(Mat);
[NL, NW] = size(Mat);
q = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isvector(x) && (min(x) > 1) && (max(mod(x,1))==0);
Chk3 = @(x) isnumeric(x) && isscalar(x) && (x > 0);
Chk4 = @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= 1);
Chkx = @(x) isnumeric(x) && ismatrix(x) && (min(size(x)) > 10);
addRequired(q,'Mat', Chkx);
addParameter(q,'m',[floor(NL/10) floor(NW/10)],Chk2);
addParameter(q,'tau',1, Chk);
addParameter(q,'r',20, Chk3);
addParameter(q,'ps',.7, Chk4);
addParameter(q,'Logx',exp(1),Chk3);
addParameter(q,'Lock',true,@(x) islogical(x));
parse(q,Mat,varargin{:})
tau = q.Results.tau;  

if (NL > 128 || NW > 128) && q.Results.Lock
    error(['To prevent memory overload, matrix width & length' ...
        ' must have <= 128 elements. \nTo estimate SampEn2D ', ...
        'for the current matrix (%dx%d) change Lock to false.', ...
        '\CAUTION: unlocking the safe matrix size may cause ', ...
        'MatLab to crash.'],NL,NW)
end
if length(q.Results.m)==1
    mL = q.Results.m; 
    mW = q.Results.m;
else
    mL = q.Results.m(1); 
    mW = q.Results.m(2);
end

NL = NL - (mL-1)*tau;
NW = NW - (mW-1)*tau;

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

Cij = -ones(p-1);
for k = 1:p-1
    Temp = abs(X(k+1:p,1:mL,1:mW) - X(k,1:mL,1:mW))<=q.Results.r;
    Cij(1:end-k+1,k) = sum(Temp,[2,3]);    
end

Dm = sum((Cij(:)/(mL*mW))>=q.Results.ps)/(p*(p-1)/2);
Esp2D = -log(Dm)/log(q.Results.Logx);
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub