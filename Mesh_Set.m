function [ grid, X ] = Mesh_Set( points, R_left, R_right, J)
%MESH_SET Initialize the mesh

h_half = (R_right - R_left) / 2 / J;
grid = ones(2,J);
grid(1,:) = R_left + h_half : h_half * 2 : R_right - h_half;
grid(2,:) = h_half * grid(2,:);
X = points * grid(2,:) + repmat(grid(1,:),size(points,1),1);
    
end

