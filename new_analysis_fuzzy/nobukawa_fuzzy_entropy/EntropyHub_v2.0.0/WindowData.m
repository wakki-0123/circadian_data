function [WinData, Log] = WindowData(Data, varargin)
% WindowData restructures a univariate/multivariate dataset into a collection of subsequence windows.
%
%   [WinData, Log] = WindowData(Data) 
% 
%   Windows the sequence(s) given in ``Data`` into a collection of subsequnces  
%   of floor(N/5) elements with no overlap, excluding any remainder
%   elements that do not fill the final window.
%   If ``Data`` is a univariate sequence (vector), ``Windata`` is a cell of 5
%   vectors. If ``Data`` is a set of multivariate sequences (NxM matrix), 
%   each of M columns is treated as a sequence with N elements 
%   and ``WinData`` is a cell of 5 matrices of size [(floor*N,5), M]. 
%   ``Log`` contains information about the windowing process, including:
%       * ``DataType``      - The type of data sequence passed as ``Data``
%       * ``DataLength``    - The number of sequence elements in ``Data``
%       * ``WindowLength``  - The number of elements in each window of ``WinData``
%       * ``WindowOverlap`` - The number of overlapping elements between windows
%       * ``TotalWindows``  - The number of windows extracted from ``Data`` 
%       * ``Mode``          - Decision to include or exclude any remaining sequence elements (``< WinLen``) that do not fill the window.
% 
%   [WinData, Log] = WindowData(Data, name, value, ...)
% 
%   Windows the sequence(s) given in ``Data`` into a collection of subsequnces
%   using the specified name/value pair arguments:
% 
%       :WinLen:  - Number of elements in each window, a positive integer (>10)
%       :Overlap: - Number of overlapping elements between windows, a positive integer (< WinLen)
%       :Mode:    - Decision to include or exclude any remaining sequence elements (< ``WinLen``) that do not fill the window, a string - either ``"include"`` or ``"exclude"`` (default).
%
%   See also:
%       ExampleData
% 

narginchk(1,7)
Data = squeeze(Data);
if isvector(Data)
    if isrow(Data); Data = Data'; end
    Log.DataType = "single univariate vector (1 sequence)";

elseif ismatrix(Data) 
    Dn = size(Data, 2);
    Log.DataType = sprintf("multivariate matrix (%d vectors)",Dn);

end

N = size(Data,1);
Chk = @(x) isnumeric(x) && ((isvector(x) && length(x)>10) || ...
               (ismatrix(x) && (size(x,1)>10) && size(x,2)>1));
Chk1 = @(x) isnumeric(x) && isscalar(x) && (mod(x,1)==0) && (x>10) && (x<N);
Chk2 = @(x) isnumeric(x) && isscalar(x) && (mod(x,1)==0) && ismember(x,0:N-1);

p = inputParser;
addRequired(p, 'Data', Chk);
addParameter(p,'WinLen', floor(N/5), Chk1);
addParameter(p,'Overlap', 0, Chk2);
addParameter(p,'Mode', "exlcude", @(x) (ischar(x) || isstring(x)) && any(strcmpi(string(x), ["include","exclude"])));
parse(p, Data,varargin{:})

Data = p.Results.Data;
WinLen = p.Results.WinLen;
Overlap = p.Results.Overlap;
Mode = p.Results.Mode;
assert(Overlap<WinLen, ...
    "The number of overlap samples must be less than the number of samples in the window!")

if WinLen<=10
    warning("There are <= 10 samples in the windowed signals." + ...
        "These signals are too short to use with any EntropyHub entropy estimator!")
end

M = floor((N - Overlap)/(WinLen - Overlap));
Step = (WinLen-Overlap);
q=1; Xout = cell(M,1);
for k = 1:Step:M*Step
    Xout{q} = Data(k:k+WinLen-1,:);
    q = q+1;
end

if strcmpi(Mode, "include") && (k+WinLen-1)<N;    Xout{q} = Data(k+Step:end,:); M=M+1; end

% WinLen+k < N
WinData = Xout;

Log.DataLength = N;
Log.WindowLength = WinLen;
Log.WindowOverlap = Overlap;
Log.TotalWindows = M;
Log.Mode = Mode;

end


%   Copyright 2024 Matthew W. Flood, EntropyHub
%   For Terms of Use see https://github.com/MattWillFlood/EntropyHub