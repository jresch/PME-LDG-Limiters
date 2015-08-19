function [ u_coord ] = Limiter( u_coord, loop, type_limiter )
%LIMITER Implement the limiters to the solution
switch type_limiter
    case 0
    case 1
        u_coord = limiter_zq(u_coord, loop);
    case 2
        u_coord = limiter_yy(u_coord, loop);
end

end

