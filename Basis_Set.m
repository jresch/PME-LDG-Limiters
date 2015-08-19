function [ psi, psi_z ] = Basis_Set( points, basis_order )

psi_full = zeros(size(points,1),5);
psi_z_full = zeros(size(points,1),5);
psi_full(:, 1) = ones(size(points,1),1);
psi_full(:, 2) = points;
psi_full(:, 3) = points .^ 2 * 1.5 - 0.5;
psi_full(:, 4) = points .^ 3 * 2.5 - points * 1.5;
psi_full(:, 5) = points .^ 4 * 4.375 - points .^ 2 * 3.75 + 0.375;

psi_z_full(:, 2) = ones(size(points,1),1);
psi_z_full(:, 3) = points * 3;
psi_z_full(:, 4) = points .^ 2 * 7.5 - 1.5;
psi_z_full(:, 5) = points .^ 3 * 17.5 - points * 7.5;

psi = psi_full(:, 1:basis_order+1);
psi_z = psi_z_full(:, 1:basis_order+1);
end

