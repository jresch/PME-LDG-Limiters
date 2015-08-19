function [ u_coord ] = limiter_yy( u_coord, loop )
%LIMITER_YY

global psi;
global track_pos;

% This limiter do not diminish the oscillation.
%% Positive
% maske sure all the cell averages >= 0
% Because of machine error, sometimes the cell average will be extreme
% small negative value.
u_coord(1, u_coord(1,:) < 0) = 0;
% some points < 0
u = psi * u_coord;
flag_pos = sum(u < 0, 1) > 0;
track_pos(flag_pos,loop) = loop;
alpha = (u_coord(1, flag_pos) - eps) ./ (u_coord(1, flag_pos) - min(u(:, flag_pos), [], 1));
alpha = max(alpha, 0);
% If only require the solution >= 0 , sometimes because of the machine
% error, there will be some extreme small negative points after adjustion.
% So the solution should >= eps after adjustion and it won't affect the
% accuracy of the solution. But sometimes the cell average is extreme
% small, then the corresponding alpha should be zero.
assert(sum(alpha > 1) == 0);
u_coord(2:end, flag_pos) = u_coord(2:end, flag_pos) * diag(alpha);

u = psi * u_coord;
assert(sum(sum(u<0)) == 0);
end

