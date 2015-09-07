function [ u_coord ] = Limiter( u_coord, loop, type_limiter )

global psi;
global track_mean;
global track_pos;

index_mean = floor(type_limiter/100);
index_osc = floor(mod(type_limiter,100)/10);
index_pos = mod(type_limiter,10);
if index_osc ~= 0
    u_coord = limiter_osc(u_coord, loop, type_limiter);
end
if index_pos ~= 0
    u_coord = limiter_pos(u_coord, loop, type_limiter);
else
    u = psi * u_coord;
    flag_mean = u_coord(1,:)' < 0;
    flag_pos = sum(u < 0, 1) > 0;
    track_mean(flag_mean,loop) = u_coord(1, flag_mean)';
    track_pos(flag_pos,loop) = loop;
    if index_mean ~= 0
        u_coord(1, flag_mean) = 0;
    end
end
end

