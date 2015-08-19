function [ u_coord ] = limiter_pos( u_coord, loop )

global psi;
global dx;
global mu;
global track_mean;
global track_osc;
global track_pos;
basis_order = size(psi, 2) - 1;
threshold = mu * dx^2;
u = psi * u_coord;
u_mean = u_coord(1,:)';
%% Oscillation
diff_left = u_mean - u(1,:)';
diff_right = u(end,:)' - u_mean;
[diff_left_mod, flag_left] = minmod(diff_left,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
[diff_right_mod, flag_right] = minmod(diff_right,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
flag_osc = (flag_left + flag_right) > 0;
track_osc(flag_osc,loop) = loop;
% diminish the oscillation
if basis_order == 1
    u_coord(2,flag_osc) = diff_left_mod(flag_osc)';
else
    if basis_order > 2
        u_coord(4:end,flag_osc) = 0;
    end
    u_coord(2,flag_osc) = (diff_right_mod(flag_osc) + diff_left_mod(flag_osc))' / 2;
    u_coord(3,flag_osc) = (diff_right_mod(flag_osc) - diff_left_mod(flag_osc))' / 2;
end

%% Positive
% maske sure all the cell averages >= 0
flag_mean = u_coord(1,:)' < 0;
track_mean(flag_mean,loop) = u_coord(1, flag_mean)';
u_coord(1, flag_mean) = 0;
% some points < 0
u = psi * u_coord;
flag_pos = sum(u < 0, 1) > 0;
track_pos(flag_pos,loop) = loop;
if basis_order ~= 1
    u_coord(3:end,flag_pos) = 0;
end
u = psi * u_coord;
flag_left = u(1,:) < 0;
u_coord(2,flag_left) = u_coord(1,flag_left);
flag_right = u(end,:) < 0;
u_coord(2,flag_right) = -u_coord(1,flag_right);

u = psi * u_coord;
assert(sum(sum(u<0)) == 0);
end

