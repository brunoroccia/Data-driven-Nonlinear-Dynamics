function [etilde, stilde] = phiES(e, s, edata, sdata, p)
%phiES Elementwise data-driven projection in (e,s) space.
%
%   [etilde, stilde] = phiES(e, s, edata, sdata, p)
%
% Computes the closest discrete data point (etilde,stilde) to a continuous
% point (e,s) according to a separable metric.
%
% Historical note (compatibility):
%   The original implementation used a mixed p/q form with q = p/(1-p).
%   In the Duffing example script p is set to 1, and the original code
%   explicitly suggested using only the e-distance in this case.
%   To preserve behavior and avoid numerical pathologies at p=1, we use:
%       p == 1:  dist = |e - edata|
%       p ~= 1:  dist = |e - edata|^p / p  +  |s - sdata|^q / q,
%                with q = p/(p-1) (Holder conjugate)
%
% Inputs
%   e, s   : scalars
%   edata  : vector of discrete elongations
%   sdata  : vector of discrete forces
%   p      : exponent parameter

edata = edata(:);
sdata = sdata(:);

if abs(p - 1) < 10*eps
    dist = abs(e - edata);
else
    q = p/(p-1);
    dist = (abs(e - edata).^p) ./ p + (abs(s - sdata).^q) ./ q;
end

[~, idx] = min(dist);
etilde = edata(idx);
stilde = sdata(idx);

end
