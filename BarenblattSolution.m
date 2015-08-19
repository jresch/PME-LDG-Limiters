function [ y ] = BarenblattSolution( x, t, m )

k = 1 / (m+1);
temp = 1 - k * (m-1) / 2 / m / t^(2 * k) * abs(x) .^ 2;
temp = max(temp, zeros(size(temp)));
y = t^-k * temp .^ (1/(m-1));

end

