function dtrSignal = detrendSample(signal, lambda, wlength)
% Detrend based on MP tarvainen et al., 2002
% Speed improved for large signals (more than 10000 samples) at expense of small
% errors at the window borders
%   INPUTS:
%       signal [*] - Real vector or matrix as input
%       lambda - lambda parameter of the detrend filter (see reference)
%       wlength - dimension in samples of the separation windows (order == 1 is equivalent to original paper code)      
% 
%   OUTPUTS:
%       dtrSignalSignal - detrended HRV signal
% 
%   REFERENCES:
%       [ 1 ] M. P. Tarvainen, P. O. Ranta-aho and P. A. Karjalainen, 
%       "An advanced detrending method with application to HRV analysis," 
%       in IEEE Transactions on Biomedical Engineering, vol. 49, no. 2, pp. 172-175, Feb. 2002.
%       doi: 10.1109/10.979357

% Output matrix allocation
S = size(signal);
dtrSignal = zeros(S);

% Last window size
wLast = mod(S(1), wlength);
nWindows = floor(S(1) / wlength);

% Coefficients creation
iMain = eye(wlength);
d2Main = spdiags(ones(wlength-2, 1)*[1 -2 1], 0:2, wlength-2, wlength);
coeffMain = (iMain - (iMain + lambda^2 * (d2Main'*d2Main)) \ iMain);
if(wLast > 0)
    iEnd = eye(wLast);
    d2End = spdiags(ones(wLast-2, 1)*[1 -2 1], 0:2, wLast-2, wLast);
    coeffEnd = (iEnd - (iEnd + lambda^2 * (d2End'*d2End)) \ iEnd);
end

% Application to main windows
for i = 1:wlength:(S(1) - wLast)
    dtrSignal(1*i:i-1+wlength,:) = coeffMain * signal(1*i:i-1+wlength,:);
end

% Application to last window
if(wLast > 0)
    dtrSignal(S(1)-wLast+1:end, :) = coeffEnd *  signal(S(1)-wLast+1:end, :);
end
end

