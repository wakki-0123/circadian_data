function [XFuzz, Ps1, Ps2] = XFuzzEn(Sig1, Sig2, varargin) 
% XFuzzEn  estimates the cross-fuzzy entropy between two univariate data sequences.
%
%   [XFuzz, Ps1, Ps2] = XFuzzEn(Sig1, Sig2) 
% 
%   Returns the cross-fuzzy entropy estimates (``XFuzz``) and the average fuzzy
%   distances (``m: Ps1``, ``m+1: Ps2``) for ``m`` = [1,2] estimated for the data sequences
%   contained in ``Sig1`` and ``Sig2``, using the default parameters: 
%   embedding dimension = 2, time delay = 1, fuzzy function (``Fx``) = ``'default'``,
%   fuzzy function paramters (``r``) = [0.2, 2], logarithm = natural
%
%   [XFuzz, Ps1, Ps2] = XFuzzEn(Sig1, Sig2, name, value, ...)
% 
%   Returns the cross-fuzzy entropy estimates (``XFuzz``) for dimensions = [1, ..., ``m``]
%   estimated for the data sequences in ``Sig1`` and ``Sig2`` using the specified name/value 
%   pair arguments:
% 
%      * ``m``     - Embedding Dimension, a positive integer   [default: 2]
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
%
%   For further information on the name/value paire arguments, see
%   the `EntropyHub guide <https://github.com/MattWillFlood/EntropyHub/blob/main/EntropyHub%20Guide.pdf>`_
%   
%   See also:
%         FuzzEn, XSampEn, XApEn, FuzzEn2D, XMSEn, MSEn
% 
%   References:
%     [1] Hong-Bo Xie, et al.,
%           "Cross-fuzzy entropy: A new method to test pattern synchrony of
%           bivariate time series." 
%           Information Sciences 
%           180.9 (2010): 1715-1724.
%
%     [2] Hamed Azami, et al.
%           "Fuzzy Entropy Metrics for the Analysis of Biomedical Signals: 
%           Assessment and Comparison"
%           IEEE Access
%           7 (2019): 104833-104847

narginchk(2,12)
p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 0) && (mod(x,1)==0);
Chk1 = @(x) isnumeric(x) && isvector(x) && numel(x)>=10;
Chk2 = @(x) isscalar(x) && (x > 0) && isnumeric(x);
Chk3 = {'sigmoid','modsampen','default','gudermannian', 'bell', 'triangular', 'trapezoidal1', 'trapezoidal2', 'z_shaped', 'gaussian', 'constgaussian'};
addRequired(p,'Sig1',Chk1);
addRequired(p,'Sig2',Chk1);
addParameter(p,'m',2,Chk);
addParameter(p,'tau',1,Chk);
addParameter(p,'r',[.2, 2.0],@(x) isnumeric(x) && isvector(x) && (min(x) >= 0));
addParameter(p,'Fx','default',@(x) (ischar(x) || isstring(x)) && any(validatestring(lower(x),Chk3)));
addParameter(p,'Logx',exp(1),Chk2);
parse(p,Sig1, Sig2, varargin{:})
m = p.Results.m; tau = p.Results.tau; 
r = p.Results.r; Logx = p.Results.Logx; 
if numel(r) == 2 && strcmpi(p.Results.Fx,'gudermannian')
    r = r(1);
    fprintf('Multiple values for r entered. First value used.') 
end

S1 = Sig1(:); S2 = Sig2(:);
N1 = length(S1); N2 = length(S2);

Fun = str2func(lower(p.Results.Fx));
m = m+1;
Sx1 = zeros(N1,m);
Sx2 = zeros(N2,m);
for k = 1:m
    Sx1(1:N1-(k-1)*tau,k) = S1(1 + (k-1)*tau:N1);
    Sx2(1:N2-(k-1)*tau,k) = S2(1 + (k-1)*tau:N2);
end

Ps1 = zeros(1,m);
Ps2 = zeros(1,m-1);
Ps1(1) = 1; %Does this need to be changed
for k = 2:m
    N1x = N1 - k*tau;
    N1y = N2 - k*tau;
    N2x = N1 - (k-1)*tau;
    N2y = N2 - (k-1)*tau;

    A = Sx1(1:N2x,1:k) - mean(Sx1(1:N2x,1:k),2);
    B = Sx2(1:N2y,1:k) - mean(Sx2(1:N2y,1:k),2);
    d2 = zeros(N2x,N2y);
    for p = 1:N2x
        Mu2 = max(abs(A(p,:) - B),[],2);
        d2(p,:) = Fun(Mu2,r);
    end    
    Ps1(k) = mean(mean(d2(1:N1x,1:N1y)));
    Ps2(k-1) = mean(mean(d2));
end
XFuzz = (log(Ps1(1:end-1)) - log(Ps2))/log(Logx);
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
% elseif r==0 && numel(x)==1
%     y = 0;
% elseif r == 1
%     y = exp(-(x - min(x)));
% else
%     error('When Fx = "Linear", r must be 0 or 1')
% %     y = -(x - max(x))/range(x);
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
