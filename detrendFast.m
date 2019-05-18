function dtrSignal = detrendFast(signal, lambda, order)
% Detrend based on MP tarvainen et al., 2002
% Speed improved for large signals (more than 10000 samples) at expense of small
% errors at the window borders
%   INPUTS:
%       signal [*] - Real vector or matrix as input
%       lambda - lambda parameter of the detrend filter (see reference)
%       order - number of divisions of the signal (order == 1 is equivalent to original paper code)      
% 
%   OUTPUTS:
%       dtrSignal - detrended HRV signal
% 
%   REFERENCES:
%       [ 1 ] M. P. Tarvainen, P. O. Ranta-aho and P. A. Karjalainen, 
%       "An advanced detrending method with application to HRV analysis," 
%       in IEEE Transactions on Biomedical Engineering, vol. 49, no. 2, pp. 172-175, Feb. 2002.
%       doi: 10.1109/10.979357

% Output matrix allocation
dtrSignal = zeros(size(signal));

% Dimensionality check, reduce order to comply with signal size
T = size(signal, 1);
if(mod(T,order)==0 && T>(2*T)/order)
    D1 = T/order;
else
    D1 = T;
end

% Coefficient generation
I = speye(D1);
D2 = spdiags(ones(D1-2,1)*[1 -2 1], 0:2, D1-2,D1);
dtrCoeff = (I-(I+lambda^2*(D2'*D2))\I);

% Detrend application
for j = 1:D1:T
    dtrSignal(1*j:j-1+D1, :) = dtrCoeff * signal(1*j:j-1+D1, :);
end