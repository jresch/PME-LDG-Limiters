function [ points, weights ] = Quadrature_Set( )
%QUADRATURE_SET Integration
% Gauss-Lobatto Quadrature, six points gauss numerical integration
% x1 = sqrt((7 - 2 * sqrt(7)) / 21);
% x2 = sqrt((7 + 2 * sqrt(7)) / 21);
% w1 = (14 + sqrt(7)) / 30;
% w2 = (14 - sqrt(7)) / 30;
% points = [-1;-x2;-x1;x1;x2;1];
% weights = [1/15,w2,w1,w1,w2,1/15];

x1 = sqrt(5-2*sqrt(10/7))/3;
x2 = sqrt(5+2*sqrt(10/7))/3;
points = [-1;-x2;-x1;0;x1;x2;1];
w1 = (322+13*sqrt(70))/900;
w2 = (322-13*sqrt(70))/900;
weights = [0,w2,w1,128/225,w1,w2,0];


end

