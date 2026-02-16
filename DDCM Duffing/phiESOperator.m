function [etilde, stilde] = phiESOperator(e, s, edata, sdata, p)
%phiESOperator ADM projection wrapper for multiple "elements".
%
%   [etilde, stilde] = phiESOperator(e, s, edata, sdata, p)
%
% Inputs
%   e, s  : vectors (Le x 1) of continuous elongations and forces
%   edata : cell array {edata{i}} with discrete elongations
%   sdata : cell array {sdata{i}} with discrete forces
%   p     : metric parameter forwarded to phiES
%
% Output
%   etilde, stilde : projected discrete pairs (Le x 1)

Le = numel(e);
etilde = zeros(Le,1);
stilde = zeros(Le,1);

for i = 1:Le
    [etilde(i), stilde(i)] = phiES(e(i), s(i), edata{i}, sdata{i}, p);
end

end
