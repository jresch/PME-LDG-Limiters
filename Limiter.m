function [ u_coord ] = Limiter( u_coord, loop, type_limiter )

index_osc = floor(mod(type_limiter,100)/10);
index_pos = mod(type_limiter,10);
if index_osc ~= 0
    u_coord = limiter_osc(u_coord, loop, type_limiter);
end
if index_pos ~= 0
    u_coord = limiter_pos(u_coord, loop, type_limiter);
end
end

