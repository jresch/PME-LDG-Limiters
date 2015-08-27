function [ u_coord ] = limiter_pos( u_coord, loop, type_limiter )

global psi;
global track_mean;
global track_pos;

index_pos = mod(type_limiter,10);
basis_order = size(psi, 2) - 1;

% maske sure all the cell averages >= 0
if type_limiter >= 100
    flag_mean = u_coord(1,:)' < 0;
    u_coord(1, flag_mean) = 0;
    if type_limiter < 200
        track_mean(flag_mean,loop) = u_coord(1, flag_mean)';
    end
end

switch index_pos
    case 1
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
    case 2
        u = psi * u_coord;
        flag_pos = sum(u < 0, 1) > 0;
        track_pos(flag_pos,loop) = loop;
        alpha = (u_coord(1, flag_pos) - eps) ./ (u_coord(1, flag_pos) - min(u(:, flag_pos), [], 1));
        alpha = max(alpha, 0);
        assert(sum(alpha > 1) == 0);
        % If only require the solution >= 0 , sometimes because of the machine
        % error, there will be some extreme small negative points after adjustion.
        % So the solution should >= eps after adjustion and it won't affect the
        % accuracy of the solution. But sometimes the cell average is extreme
        % small, then the corresponding alpha should be zero.
        u_coord(2:end, flag_pos) = u_coord(2:end, flag_pos) * diag(alpha);
end

u = psi * u_coord;
assert(sum(sum(u<0)) == 0);


end

