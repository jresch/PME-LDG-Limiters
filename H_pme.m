function [ flux_ur, q ] = H_pme( u, loop, psi, psi_z )

global m;
global weights;
global dx;
Au = @(u) m * u.^(m-1);
Gu = @(u) 2*sqrt(m)/(m+1) * u.^((m+1)/2);

% calculalte q
if mod(loop, 2)
    flux_qr = -Gu(u(end,:));
else
    flux_qr = -Gu(u(1,[2:end,1]));
end
flux_ql = flux_qr([end,1:end-1]);
A = psi' * diag(weights) * psi;
RHS = psi_z' * diag(weights) * -Gu(u) ...
    - psi(end,:)' * flux_qr ...
    + psi(1,:)' * flux_ql;
q_coord = A \ (RHS / (dx/2));
q = psi * q_coord;
% calculate u
flux_ur = ( Gu(u(1,[2:end,1])) - Gu(u(end,:)) ) ./ ( u(1,[2:end,1]) - u(end,:));
index_nan = isnan(flux_ur);
flux_ur(index_nan) = sqrt(Au(u(end,index_nan)));
if mod(loop, 2)
    flux_ur = -flux_ur .* q(1,[2:end,1]);
else
    flux_ur = -flux_ur .* q(end,:);
end
end

