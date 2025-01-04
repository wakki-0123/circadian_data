function [Fuzz2D] = FuzzEn2D(Mat, varargin)
% FuzzEn2D  estimates the bidimensional fuzzy entropy of a data matrix.
%
%   [Fuzz2D] = FuzzEn2D(Mat) 
% 
%   Returns the bidimensional fuzzy entropy estimate (``Fuzz2D``) estimated for 
%   the data matrix (``Mat``) using the default parameters: time delay = 1,
%   fuzzy function (``Fx``) = ``'default'``, fuzzy function parameters (``r``) = [0.2,2],
%   logarithm = natural, template matrix size = [floor(H/10) floor(W/10)]  
%   (where H and W represent the height (rows) and width (columns) of the data matrix ``Mat``) 
%   ** The minimum number of rows and columns of ``Mat`` must be > 10.
%
%   [Fuzz2D] = FuzzEn2D(Mat, name, value, ...)
% 
%   Returns the bidimensional fuzzy entropy (``Fuzz2D``) estimates for the data
%   matrix (``Mat``) using the specified name/value pair arguments:
% 
%      * ``m``     - Template submatrix dimensions, an integer scalar (for sub-
%        matrix with same height and width) or a two-element vector of
%        integers [height, width] with a minimum value > 1.
%        (default: [floor(H/10) floor(W/10)])
%      * ``tau``   - Time Delay, a positive integer        [default: 1]
%      * ``Fx``    - Fuzzy function name, one of the following strings:
%        {``'sigmoid'``, ``'modsampen'``, ``'default'``, ``'gudermannian'``, ``'bell'``, ``'triangular'``, ``'trapezoidal1'``, ``'trapezoidal2'``, ``'z_shaped'``, ``'gaussian'``, ``'constgaussian'``}
%      * ``r``     - Fuzzy function parameters, a 1 element scalar or a 2 element
%        vector of positive values. The ``r`` parameters for each fuzzy
%        function are defined as follows:  (default: [.2 2])
%                 sigmoid:     
%                             * r(1) = divisor of the exponential argument
%                             * r(2) = value subtracted from argument (pre-division)
%                 modsampen:    
%                             *  r(1) = divisor of the exponential argument
%                             *  r(2) = value subtracted from argument (pre-division)
%                 default:      
%                             *  r(1) = divisor of the exponential argument
%                             *  r(2) = argument exponent (pre-division)
%                 gudermannian: 
%                             *  r  = a scalar whose value is the numerator of
%                               argument to gudermannian function:
%                               GD(x) = atan(tanh(r/x)). GD(x) is 
%                               normalised to have a maximum value of 1.
%                 [DEPRICATED] linear:       
%                              * r  = an integer value. When ``r == 0``, the
%                               argument of the exponential function is 
%                               normalised between [0 1]. When ``r == 1``,
%                               the minimuum value of the exponential 
%                               argument is set to 0.                       
%                 triangular:  
%                             * r = a positive scalar whose value is the threshold 
%                               (corner point) of the triangular function.
%                 trapezoidal1:  
%                             * r = a positive scalar whose value corresponds
%                               to the upper (2r) and lower (r) corner points of the trapezoid.
%                 trapezoidal2:  
%                             * r(1) = a value corresponding to the upper corner point of the trapezoid.
%                             * r(2) = a value corresponding to the lower corner point of the trapezoid.
%                 z_shaped:  
%                             * r = a scalar whose value corresponds to the
%                               upper (2r) and lower (r) corner points of the z-shape.
%                 bell:  
%                             * r(1) = divisor of the distance value
%                             * r(2) = exponent of generalized bell-shaped function
%                 gaussian:  
%                             * r = a scalar whose value scales the slope of the Gaussian curve.
%                 constgaussian:  
%                             * r = a scalar whose value defines the lower 
%                               threshod and shape of the Gaussian curve.                         
%      * ``Logx``  - Logarithm base, a positive scalar  [default: natural]
%      * ``Lock``  - By default, ``FuzzEn2D`` only permits matrices with a maximum
%        size of 128 x 128 to prevent RAM overload. 
%        e.g. For ``Mat`` = [200 x 200], ``m = 3``, and ``tau = 1``, ``FuzzEn2D`` 
%        creates a vector of 753049836 elements. To enable matrices
%        greater than [128 x 128] elements, set ``'Lock' = false``. (default: true)
% 
%        **WARNING: unlocking the permitted matrix size may cause memory
%        errors that could lead Matlab to crash.**
%
%   See also:
%       SampEn2D, DistEn2D, FuzzEn, XFuzzEn, XMSEn
% 
%   References:
%     [1] Luiz Fernando Segato Dos Santos, et al.,
%           "Multidimensional and fuzzy sample entropy (SampEnMF) for
%           quantifying H&E histological images of colorectal cancer."
%           Computers in biology and medicine 
%           103 (2018): 148-160.
% 
%     [2] Mirvana Hilal and Anne Humeau-Heurtier,
%           "Bidimensional fuzzy entropy: Principle analysis and biomedical
%           applications."
%           41st Annual International Conference of the IEEE (EMBC) Society
%           2019.
% 
%     [3] Hamed Azami, et al.
%           "Fuzzy Entropy Metrics for the Analysis of Biomedical Signals: 
%           Assessment and Comparison"
%           IEEE Access
%           7 (2019): 104833-104847

narginchk(1,13)
Mat = squeeze(Mat);
[NL, NW] = size(Mat);
q = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk2 = @(x) isvector(x) && (min(x) > 1) && (max(mod(x,1))==0);
Chk3 = @(x) isnumeric(x) && isvector(x) && (x(1) >= 0);
Chkx = @(x) isnumeric(x) && ismatrix(x) && (min(size(x)) > 10);
Chky = {'sigmoid','modsampen','default','gudermannian', 'bell', 'triangular', 'trapezoidal1', 'trapezoidal2', 'z_shaped', 'gaussian', 'constgaussian'};
addRequired(q,'Mat', Chkx);
addParameter(q,'m',[floor(NL/10) floor(NW/10)],Chk2);
addParameter(q,'tau',1, Chk);
addParameter(q,'Fx','default',@(x) (ischar(x) || isstring(x)) && any(validatestring(lower(x),Chky)));
addParameter(q,'r', [.2*std(Mat(:),1) 2.0], Chk3);
addParameter(q,'Logx',exp(1),Chk3);
addParameter(q,'Lock',true,@(x) islogical(x));
parse(q,Mat,varargin{:})
tau = q.Results.tau; r = q.Results.r;

if numel(r) == 2 && strcmpi(q.Results.Fx,'gudermannian')
    r = r(1);
    fprintf('Multiple values for r entered. First value used.') 
end

if (NL > 128 || NW > 128) && q.Results.Lock
        error(['To prevent memory errors, matrix width & length' ...
            ' must have <= 128 elements. \nTo estimate FuzzEn2D ', ...
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
Fun = str2func(lower(q.Results.Fx));
NL = NL - mL*tau;
NW = NW - mW*tau;
X1 = zeros(NL*NW,mL,mW);
X2 = zeros(NL*NW,mL+1,mW+1);
p = 0;
for k = 1:NL
    for n = 1:NW
        p = p+1;
        Temp2 = Mat(k:tau:(mL)*tau+k,n:tau:(mW)*tau+n);
        Temp1 = Temp2(1:end-1,1:end-1);
        X1(p,:,:) = Temp1 - mean(mean(Temp1));
        X2(p,:,:) = Temp2 - mean(mean(Temp2));
    end
end

if p ~= NL*NW
    warning('Potential error with submatrix division.')
end
Ny = p*(p-1)/2;
if Ny > 300000000
    warning('Number of pairwise distance calculations is %d', Ny)
end

Y1 = zeros(1,p-1);
Y2 = zeros(1,p-1);
for k = 1:p-1
    Temp1 = max(abs(X1(k+1:end,:,:) - X1(k,:,:)),[],[2,3]);
    Y1(k) = sum(Fun(Temp1, r));
    Temp2 = max(abs(X2(k+1:end,:,:) - X2(k,:,:)),[],[2,3]);
    Y2(k) = sum(Fun(Temp2, r));
end
Fuzz2D = -log(sum(Y2)/sum(Y1))/log(q.Results.Logx);
end

function [y] = sigmoid(x,r)
if numel(r) == 1
    error('When Fx = "Sigmoid", r must be a two-element vector.')
end
y = 1./(1+exp((x-r(2))/r(1)));
end
function [y] = modsampen(x,r)
if numel(r) == 1
    error('When Fx = "Modsampen", r must be a two-element vector.')
end
y = 1./(1+exp((x-r(2))/r(1)));
end
function [y] = default(x,r)
if numel(r) == 1
    error('When Fx = "Default", r must be a two-element vector.')
end
y = exp(-(x.^r(2))/r(1));
end
function [y] = gudermannian(x,r)
y = atan(tanh(r(1)./x));
y = y/max(y);
end

% function [y] = linear(x,r)
% if r == 0 && numel(x)>1
%     y = exp(-(x - min(x))/range(x));
% elseif r == 0 && numel(x) == 1
%     y = 0;
% elseif r == 1
%     y = exp(-(x - min(x)));
% else
%     error('When Fx = "Linear", r must be 0 or 1')
% end
% end

function [y] = triangular(x,r)
    if numel(r) > 1
        error('When Fx = "Triangular", r must be a scalar > 0.')
    end
    y = 1 - (x/r);
    y(x > r) = 0;
end
function [y] = trapezoidal1(x,r)
    if numel(r) > 1
        error('When Fx = "Trapezoidal1", r must be a scalar > 0.')
    end
    y = zeros(size(x));
    y(x <= r*2) = 2 - (x(x <= r*2)/r);
    y(x <= r) = 1;
end
function [y] = trapezoidal2(x,r)
    if numel(r) ~=2 
        error('When Fx = "Trapezoidal2", r must be a two-element vector.')
    end
    y = zeros(size(x));
    y(x <= max(r)) = 1 - (x(x <= max(r))-min(r))/range(r);
    y(x <= min(r)) = 1;
end
function [y] = z_shaped(x, r)
    if numel(r) > 1
        error('When Fx = "Z_shaped", r must be a scalar > 0.')
    end
    y = zeros(size(x));

    y(x <= 2*r) = 2*(((x(x <= 2*r) - 2*r)/r).^2);
    y(x <= 1.5*r) = 1 - (2*(((x(x <= 1.5*r) - r)/r).^2));
    y(x <= r) = 1;
end
function [y] = bell(x, r)
    if numel(r) ~=2 
        error('When Fx = "Bell", r must be a two-element vector.')
    end
    y = 1./(1 + abs(x/r(1)).^(2*r(2)));
end
function [y] = gaussian(x, r)
    if numel(r) > 1
        error('When Fx = "Gaussian", r must be a scalar > 0.')
    end
    y = exp(-((x.^2)/(2*(r^2))));
end
function [y] = constgaussian(x, r)
    if numel(r) > 1
        error('When Fx = "ConstGaussian", r must be a scalar > 0.')
    end
    y = ones(size(x));
    y(x > r) = exp(-log(2)*((x(x > r) - r)/r).^2);
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub
