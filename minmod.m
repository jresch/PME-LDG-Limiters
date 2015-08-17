function [ a, flag ] = minmod( a1, a2, a3, threshold )
%MINMOD 

a = zeros(size(a1));
flag = zeros(size(a1));
% Situation 0
index0 = abs(a1) <= threshold;
a(index0) = a1(index0);
% Situation 1
index1 = (sign(a1)==sign(a2)) .* (sign(a1)==sign(a3)) .* (index0 == 0);
index1 = index1 > 0;
temp = sign(a1) .* min(abs([a1,a2,a3]),[],2);
a(index1) = temp(index1);
flag(index1) = 1;
% Situation 2
index2 = (1 - index0 - index1)>0;
flag(index2) = 2;

end

