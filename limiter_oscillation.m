function [ u_coord ] = limiter_oscillation( u_coord, loop )
%LIMITER_OSCILLATION
global psi;
global grid;
global mu;
global count_osc;

basis_order = size(psi, 2) - 1;
threshold = mu * grid(2,1)^2 * 4;
u = psi * u_coord;
u_mean = u_coord(1,:)';
diff_left = u_mean - u(1,:)';
diff_right = u(end,:)' - u_mean;
[diff_left_mod, flag_left] = minmod(diff_left,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
[diff_right_mod, flag_right] = minmod(diff_right,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
u_left_mod = u_mean - diff_left_mod;
u_right_mod = diff_right_mod + u_mean;
flag_mod = (flag_left + flag_right) > 0;
count_osc(flag_mod,loop) = loop;
if basis_order == 1
    u_coord(2,flag_mod) = diff_left_mod(flag_mod);
else
    A = psi([1,end],2:3);
    RHS = [u_left_mod(flag_mod)';u_right_mod(flag_mod)'] - psi([1,end],1) * u_coord(1,flag_mod);
    u_coord(:,flag_mod) = 0;
    u_coord(1,flag_mod) = u_mean(flag_mod)';
    u_coord(2:3,flag_mod) = A \ RHS;
end


end

