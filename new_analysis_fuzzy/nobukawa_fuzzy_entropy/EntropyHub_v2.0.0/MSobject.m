function [Mobj] = MSobject(EnType,varargin)
% MSobject  creates an object to store multiscale entropy parameters.
%
%   [Mobj] = MSobject() 
% 
%   Returns a multiscale entropy object (``Mobj``) based on that orignially
%   proposed by Costa et al. (2002) using the following default
%   parameters:
%   EnType = ``'SampEn'``, embedding dimension = 2, time delay = 1, 
%   radius = 0.2*SD(``Sig``), logarithm = natural
%    
%   [Mobj] = MSobject(EnType)
% 
%   Returns a multiscale entropy object using the specified entropy method
%   (``EnType``) and the default parameters for that entropy method.
%   To see the default parameters for a particular entropy method,          
%   type:   help `EnType`   (e.g.  ``help SampEn``)
%
%   [Mobj] = MSobject(EnType, name, value, ...)
% 
%   Returns a multiscale entropy object using the specified entropy method
%   (``EnType``) and the name/value parameters for that particular method.
%   To see the default parameters for a particular entropy method,          
%   type:   help `EnType`   (e.g.    ``help SampEn``)
%
%   ``EnType`` can be any of the following (case sensitive) entropies:
% 
%   :Base Entropies:
%   ---------------
% 
%      * ``'ApEn'``     - Approximate Entropy
%      * ``'SampEn'``    - Sample Entropy
%      * ``'FuzzEn'``    - Fuzzy Entropy
%      * ``'K2En'``      - Kolmogorov Entropy
%      * ``'PermEn'``    - Permutation Entropy	
%      * ``'CondEn'``    - Conditional Entropy	
%      * ``'DistEn'``    - Distribution Entropy	
%      * ``'DispEn'``    - Dispersion Entropy	
%      * ``'SpecEn'``    - Spectral Entropy
%      * ``'SyDyEn'``    - Symbolic Dynamic Entropy	
%      * ``'IncrEn'``    - Increment Entropy	
%      * ``'CoSiEn'``    - Cosine Similarity Entropy	
%      * ``'PhasEn'``    - Phase Entropy	
%      * ``'SlopEn'``    - Slope Entropy
%      * ``'BubbEn'``    - Bubble Entropy	
%      * ``'GridEn'``    - Grid Distribution Entropy	
%      * ``'EnofEn'``    - Entropy of Entropy	
%      * ``'AttnEn'``    - Attention Entropy
%      * ``'DivEn'``     - Diversity Entropy
%      * ``'RangEn'``    - Range Entropy
%   
%   :Cross Entropies:
%   ----------------
% 
%      * ``'XApEn'``     - Cross-Approximate Entropy
%      * ``'XSampEn'``   - Cross-Sample Entropy
%      * ``'XFuzzEn'``   - Cross-Fuzzy Entropy
%      * ``'XK2En'``     - Cross-Kolmogorov Entropy
%      * ``'XPermEn'``   - Cross-Permutation Entropy
%      * ``'XCondEn'``   - Cross-Conditional Entropy
%      * ``'XDistEn'``   - Cross-Distribution Entropy
%      * ``'XSpecEn'``   - Cross-Spectral Entropy
%
%   :Multivariate Entropies:
%   ------------------------
% 
%      * ``'MvSampEn2D'``   - Multivariate Sample Entropy
%      * ``'MvFuzzEn2D'``   - Multivariate Fuzzy Entropy
%      * ``'MvDispEn2D'``   - Multivariate Dispersion Entropy
%      * ``'MvCoSiEn2D'``   - Multivariate Cosine Similarity Entropy
%      * ``'MvPermEn2D'``   - Multivariate Permutation Entropy
% 
%   See also:
%       MSEn, XMSEn, MvMSEn, cMSEn, rMSEn, hMSEn,rXMSEn, cXMSEn, hXMSEn
%   


narginchk(0,19)
if nargin == 0
    EnType = 'SampEn';
end
p = inputParser;
Chk = {'ApEn';'SampEn';'FuzzEn';'K2En';'PermEn';'CondEn';'DistEn'; ...
    'DispEn';'SyDyEn';'IncrEn';'CoSiEn';'PhasEn';'SpecEn';'SlopEn'; ...
    'GridEn';'BubbEn';'EnofEn';'AttnEn';'DivEn';'RangEn';'XApEn';...
    'XSampEn';'XFuzzEn';'XPermEn';'XCondEn';'XDistEn';'XSpecEn';'XK2En';...
    'MvSampEn';'MvFuzzEn';'MvCoSiEn';'MvPermEn';'MvDispEn'};
addRequired(p,'EnType',@(x) any(strcmp(string(x),Chk)));
parse(p,EnType); 
assert(mod(length(varargin),2)==0)
Mobj = struct('Func',str2func(p.Results.EnType));
for k = 1:2:length(varargin)
    Mobj.(varargin{k}) = varargin{k+1};
end

Fields = fieldnames(Mobj);
if length(Fields)>1
cx =[find(strcmp(Fields, 'm')); find(strcmp(Fields, 'tau')); ...
    find(strcmp(Fields, 'r')); find(strcmp(Fields, 'c'))];
cy = [1:length(Fields)];  cy([1 cx'])=[]; 
cy = [1 cx' cy];
Mobj = orderfields(Mobj,cy);
end
end

%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub