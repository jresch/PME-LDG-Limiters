function [ u_coord ] = limiter_positive( u_coord, loop )
%LIMITER_POSITIVE 
global psi;
global grid;
global mu;
global count_pos_mean;
global count_pos_side;
global count_pos_diff;

basis_order = size(psi, 2) - 1;
threshold = mu * grid(2,1)^2 * 4;
u = psi * u_coord;
% u_mean >= 0
flag_mean = u_coord(1,:)' < 0;
count_pos_mean(flag_mean,loop) = u_coord(1, flag_mean)';
u_coord(1, flag_mean) = 0;
u_mean = u_coord(1,:)';
% diff
diff_left = u_mean - u(1,:)';
diff_right = u(end,:)' - u_mean;
[diff_left_mod, flag_left] = minmod(diff_left,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
[diff_right_mod, flag_right] = minmod(diff_right,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
flag_side = (flag_left + flag_right) > 0;
count_pos_side(flag_side,loop) = loop;
% positive
flag_left = diff_left_mod > u_mean;
diff_left_mod(flag_left) = u_mean(flag_left);
flag_right = diff_right_mod < -u_mean;
diff_right_mod(flag_right) = -u_mean(flag_right);
flag_mod = flag_side | flag_left | flag_right;
if basis_order == 1
    u_coord(2,flag_mod) = min(abs(diff_left_mod(flag_mod)), abs(diff_right_mod(flag_mod))) .* sign(diff_left_mod(flag_mod));
else
    temp1 = diff_right_mod + diff_left_mod;
    temp2 = diff_right_mod - diff_left_mod;
    x_extreme = - temp1 / 3 ./ temp2;
    flag_diff = (x_extreme < 1) & (x_extreme > -1) & (diff_left_mod < diff_right_mod);
    alpha = 3 * temp2 .* u_mean ./ (temp1 .^ 2 / 4 + 3/4 * temp2 .^ 2);
    flag_diff = flag_diff & (alpha < 1);
    count_pos_diff(flag_diff | flag_left | flag_right,loop) = loop;
     
    diff_left_mod(flag_diff) = diff_left_mod(flag_diff) .* alpha(flag_diff);
    diff_right_mod(flag_diff) = diff_right_mod(flag_diff) .* alpha(flag_diff);
    flag_mod = flag_mod | flag_diff;
    u_coord(2:end,flag_mod) = 0;
    u_coord(2,flag_mod) = (diff_right_mod(flag_mod)' + diff_left_mod(flag_mod)') / 2;
    u_coord(3,flag_mod) = (diff_right_mod(flag_mod)' - diff_left_mod(flag_mod)') / 2;
end
% should try to limite the order of solution when there is negative points
% calculate the negative mean, try to prove that in PME, it is useless to
% adjust the mean value, (theory proof)

end