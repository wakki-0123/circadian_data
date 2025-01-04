function [XSpec, BandEn] = XSpecEn(Sig1, Sig2, varargin)
% XSpecEn  estimates the cross-spectral entropy between two univariate  data sequences.
%
%   [XSpec, BandEn] = XSpecEn(Sig1, Sig2) 
% 
%   Returns the cross-spectral entropy estimate (``XSpec``) of the full cross-
%   spectrum and the within-band entropy (``BandEn``) estimated between the data 
%   sequences contained in ``Sig1`` and ``Sig2`` using the default  parameters: 
%   N-point FFT = 2 * max length of ``Sig1``/``Sig2``,
%   normalised band edge frequencies = [0 1], logarithm = natural, 
%   normalisation = w.r.t # of spectrum/band frequency values.
%
%   [XSpec, BandEn] = XSpecEn(Sig1, Sig2, name, value, ...)
% 
%   Returns the cross-spectral entropy (``XSpec``) and the within-band entropy 
%   (``BandEn``) estimate between the data sequences contained in  ``Sig1`` 
%   and ``Sig2`` using the following specified name/value pair arguments:
% 
%       * ``N``     - Resolution of spectrum (N-point FFT), an integer > 1
%       * ``Freqs`` - Normalised band edge frequencies, a scalar in range [0 1]
%         where 1 corresponds to the Nyquist frequency (Fs/2).
%         *Note: When no band frequencies are entered, ``BandEn == SpecEn``
%       * ``Logx``  - Logarithm base, a positive scalar     [default: natural]
%       * ``Norm``  - Normalisation of ``XSpec`` value, a boolean:
%              *   [false]  no normalisation.
%              *   [true]   normalises w.r.t # of frequency values within the 
%                  spectrum/band   [default]
%
%   See also:
%         SpecEn, fft, XDistEn, periodogram, XSampEn, XApEn
%  
%   References:
%     [1]   Matthew W. Flood,
%               "XSpecEn - EntropyHub Project"
%               (2021) https://github.com/MattWillFlood/EntropyHub
% 

narginchk(2,10)
p = inputParser;
Chk = @(x) isnumeric(x) && isscalar(x) && (x > 1) && (mod(x,1)==0);
Chk1 = @(x) isnumeric(x) && isvector(x) && numel(x)>=10;
Chk2 = @(x) isvector(x) && (min(x) >= 0) && (max(x) <= 1) && (length(x)==2);
addRequired(p,'Sig1',Chk1);
addRequired(p,'Sig2',Chk1);
addParameter(p,'N',2*max(numel(Sig1),numel(Sig2)) + 1,Chk);
addParameter(p,'Freqs',[0 1],Chk2);
addParameter(p,'Logx',exp(1),@(x) isscalar(x) && (x >= 0));
addParameter(p,'Norm',true,@(x) islogical(x));
parse(p,Sig1, Sig2,varargin{:})
N = p.Results.N; Freqs = p.Results.Freqs; 
Logx = p.Results.Logx; Norm = p.Results.Norm;

S1 = Sig1(:); S2 = Sig2(:);
if Logx == 0,  Logx = exp(1);  end
% if N == 0 || isempty(N),  N = length(S1); end
Fx = ceil(N/2);
Freqs = round(Freqs*Fx);
Freqs(Freqs==0) = 1;

if length(Freqs) < 2
    error('Freqs must contain at least two values.')
elseif Freqs(1) > Freqs(2)
    error('Lower band frequency must come first.')
elseif diff(Freqs) < 1
    error('Spectrum resoution too low to determine bandwidth.') 
elseif min(Freqs)<0 || max(Freqs)> Fx
    error('Freqs must be normalized w.r.t sampling frequency [0 1].')    
end

Pt = abs(fft(conv(S1,S2),N));
Pxx = Pt(1:Fx)/sum(Pt(1:Fx));
XSpec = -(Pxx(:)'*log(Pxx(:)))/log(Logx);
Pband = (Pxx(Freqs(1):Freqs(2)))/sum(Pxx(Freqs(1):Freqs(2)));
BandEn = -(Pband*log(Pband)')/log(Logx);

if Norm
    XSpec = XSpec/(log(Fx)/log(Logx));
    BandEn = BandEn/(log(diff(Freqs)+1)/log(Logx));
end
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub