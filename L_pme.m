function [ ut_coord ] = L_pme( u_coord, loop, type_limiter )

global m;
global c;
global p;
global weights;
global dx;
global dt;
global psi;
global psi_z;
global track_mean;

Au = @(u) m * u.^(m-1);
u = psi * u_coord;
assert(sum(isnan(u(:))) == 0);
if type_limiter == 2
    [flux_ur_linear,~] = H_pme(u, loop, psi(:,1), psi_z(:,1));
    [flux_ur,q] = H_pme(u, loop, psi, psi_z);

    lambda = dt / dx;
    gamma = eps - u_coord(1,:) ...
        + lambda * (flux_ur_linear - flux_ur_linear([end,1:end-1]));
    gamma = gamma / lambda;
    F_r = flux_ur - flux_ur_linear;
    F_l = F_r([end,1:end-1]);
    theta_r = ones(size(F_r));
    theta_l = theta_r;
    % case0
    unew_coord = u_coord(1,:) - lambda * (flux_ur - flux_ur([end,1:end-1]));
    case0 = unew_coord < 0;
    track_mean(case0,loop) = unew_coord(case0)';
    % case1 F_r <= 0, F_l >= 0
    % case2 F_r <= 0, F_l < 0
    case2 = (F_r <= 0) & (F_l < 0) & case0;
    theta_l(case2) = min(1, (F_r(case2) + gamma(case2)) ./ F_l(case2));
    % case3 F_r > 0, F_l >= 0
    case3 = (F_r > 0) & (F_l >= 0) & case0;
    theta_r(case3) = min(1, (F_l(case3) - gamma(case3)) ./ F_r(case3));
    % case4 F_r > 0, F_l < 0
    case4 = (F_r > 0) & (F_l < 0) & case0;
    theta_r(case4) = min(1, - gamma(case4) ./ (F_r(case4) - F_l(case4)));
    theta_l(case4) = theta_r(case4);
    
    theta_r = min(theta_r, theta_l([2:end,1]));
    flux_ur = theta_r .* F_r + flux_ur_linear;
else
    [flux_ur,q] = H_pme(u, loop, psi, psi_z);
end
flux_ul = flux_ur([end,1:end-1]);

RHS = psi_z' * diag(weights) * (-sqrt(Au(u)) .* q) ...
    - psi(end,:)' * flux_ur ...
    + psi(1,:)' * flux_ul ...
    - dx/2 * c * psi' * diag(weights) * (u .^ p);
A = psi' * diag(weights) * psi;
ut_coord = A \ (RHS / (dx/2));

unew_coord = u_coord + dt * ut_coord;

assert(sum(unew_coord(1,:) < -eps*10) == 0);

end