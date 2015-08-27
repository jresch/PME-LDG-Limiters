function [ u_coord ] = limiter_osc( u_coord, loop, type_limiter )

global psi;
global dx;
global mu;
global weights;
global track_osc;

basis_order = size(psi, 2) - 1;
index_osc = floor(mod(type_limiter,100)/10);
threshold = mu * dx^2;
u = psi * u_coord;
u_mean = u_coord(1,:)';
% discover the oscillation
diff_left = u_mean - u(1,:)';
diff_right = u(end,:)' - u_mean;
[diff_left_mod, flag_left] = minmod(diff_left,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
[diff_right_mod, flag_right] = minmod(diff_right,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
flag_osc = (flag_left + flag_right) > 0;
track_osc(flag_osc,loop) = loop;
% diminish the oscillation
switch index_osc
    case 1
        if basis_order > 1
            u_coord(3:end,flag_osc) = 0;
            u = psi * u_coord;
            u_mean = u_coord(1,:)';
            diff_left = u_mean - u(1,:)';
            diff_right = u(end,:)' - u_mean;
            [diff_left_mod, flag_left] = minmod(diff_left,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
            % [diff_right_mod, flag_right] = minmod(diff_right,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
            flag_osc = (flag_left > 0) & flag_osc;
            track_osc(flag_osc,loop) = -loop;
        end
        u_coord(2,flag_osc) = diff_left_mod(flag_osc)';
    case 2
        if basis_order == 1
            u_coord(2,flag_osc) = diff_left_mod(flag_osc)';
        else
            if basis_order > 2
                u_coord(4:end,flag_osc) = 0;
                u = psi * u_coord;
                u_mean = u_coord(1,:)';
                diff_left = u_mean - u(1,:)';
                diff_right = u(end,:)' - u_mean;
                [diff_left_mod, flag_left] = minmod(diff_left,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
                [diff_right_mod, flag_right] = minmod(diff_right,u_mean([2:end,1])-u_mean,u_mean-u_mean([end,1:end-1]),threshold);
                flag_osc = ((flag_left + flag_right) > 0) & flag_osc;
                track_osc(flag_osc,loop) = -loop;
            end
            u_coord(2,flag_osc) = (diff_right_mod(flag_osc) + diff_left_mod(flag_osc))' / 2;
            u_coord(3,flag_osc) = (diff_right_mod(flag_osc) - diff_left_mod(flag_osc))' / 2;
        end
    case 3
        W = 1 ./ sqrt(diag(psi' * diag(weights) * psi));
        W = diag(W(2:end));
        A = [psi(end,2:end);psi(1,2:end)];
        B = [diff_right_mod(flag_osc)';-diff_left_mod(flag_osc)'] - A * u_coord(2:end, flag_osc);
        u_coord(2:end,flag_osc) = W * ((A * W) \ B) + u_coord(2:end, flag_osc);
end

end

