function [ err ] = PME( PARA, uI, path_data, path_report )

global m;
global c;
global p;
global weights;
global dx;
global dt;
global psi;
global psi_z;
global mu;
global track_mean;
global track_osc;
global track_pos;

% extract parameters
m = PARA(1,1);
c = PARA(1,2);
p = PARA(1,3);
R_left = PARA(2,1);
R_right = PARA(2,2);
dT = PARA(2,3);
J = PARA(3,1);
basis_order = PARA(3,2);
type_problem = PARA(4,1);
type_limiter = PARA(4,2);
% mu = 4/3*(m-2)/m/(m^2-1)/T0^(3/(m+1));
mu = 1;
num_fig = 200;
% Quadrature
[points,weights] = Quadrature_Set();
% Mesh
[grid,X] = Mesh_Set(points,R_left,R_right,J);
% Space & Time Step
dx = max(grid(2,:)) * 2;
if type_problem == 0
    CFL =  PARA(4,3);
else
    CFL_list = [ 50, 70, 110, 130;
            170, 270, 390, 490;
            490, 790, 1110, 1450;
            1110, 1850, 2610, 3370;];
    CFL = CFL_list(basis_order,max(1,floor(m/2)));
end
loops = ceil(dT / dx / dx * CFL);
dt = dT / loops;
% Basis
[psi,psi_z] = Basis_Set(points,basis_order);
% Limiters
track_mean = sparse(J,loops);
track_osc = sparse(J,loops);
track_pos = sparse(J,loops);
% initialize
u = uI(X);
A = psi' * diag(weights) * psi;
RHS = psi' * diag(weights) * u;
u_coord = A \ RHS;
u_coord = Limiter(u_coord,1,type_limiter);
% Solve
for l = 1:loops
    u0_coord = u_coord;
    % step 1
    u_coord = L_pme(u0_coord,l,type_limiter);
    u1_coord = u0_coord + u_coord * dt;
    u1_coord = Limiter(u1_coord,l,type_limiter);
    % step 2
    u_coord = L_pme(u1_coord,l,type_limiter);
    u2_coord = u0_coord * 3/4 + u1_coord /4 + u_coord * dt/4;
    u2_coord = Limiter(u2_coord,l,type_limiter);
    % step 3
    u_coord = L_pme(u2_coord,l,type_limiter);
    u3_coord = u0_coord / 3 + u2_coord * 2/3 + u_coord * dt*2/3;
    u3_coord = Limiter(u3_coord,l,type_limiter);
    u_coord = u3_coord;
    if (type_problem == 2) && (mod(l,floor(loops/num_fig)) == 1)
        u = psi * u_coord;
        filename = sprintf('PME-m%.1f-c%.1f-p%.1f-basis%d-lim%d-J%d-N%d',m,c,p,basis_order,type_limiter,J,ceil(l/floor(loops/num_fig)));
        figname = sprintf('%s%s.png',path_data,filename);
        filename = sprintf('%s%s.mat',path_data,filename);
        save(filename,'u','u_coord','X','R_left','R_right','l','loops');
        fig = figure;
        set(fig,'visible','off');
        plot(X(:),u(:),'.-');
        titlename = sprintf('m = %d, J = %d, order = %d, lim = %.3d, N = %d', m, J, basis_order, type_limiter,ceil(l/floor(loops/num_fig)));
        title(titlename);
        axis([R_left,R_right,-0.1,1.2]);
        print(fig,'-dpng',figname);
    end
end

if type_problem <= 1
    u = psi * u_coord;
    u_exact_pme = BarenblattSolution(X, 1+dT, m);
    err = sqrt(dx/2 * sum(weights * (u_exact_pme - psi * u_coord) .^ 2));
    if ~isnan(err)
        filename = sprintf('PME-m%.1f-c%.1f-p%.1f-basis%d-lim%.3d-J%d',m,c,p,basis_order,type_limiter,J);
        figname = sprintf('%s%s.png',path_data,filename);
        filename = sprintf('%s%s.mat',path_data,filename);
        save(filename,'u','u_coord','u_exact_pme','X','R_left','R_right','loops');
        fig = figure;
        set(fig,'visible','off');
        plot(X(:),u(:),'b.-',X(:),u_exact_pme(:),'r');
        titlename = sprintf('m = %d, J = %d, order = %d, lim = %.3d', m, J, basis_order, type_limiter);
        title(titlename);
        axis([R_left,R_right,-0.1,1.2]);
        print(fig,'-dpng',figname);
    end
elseif type_problem == 2
    u = psi * u_coord;
    filename = sprintf('PME-m%.1f-c%.1f-p%.1f-basis%d-lim%.3d-J%d-N%d',m,c,p,basis_order,type_limiter,J,ceil(loops/floor(loops/num_fig))+1);
    figname = sprintf('%s%s.png',path_data,filename);
    filename = sprintf('%s%s.mat',path_data,filename);
    save(filename,'u','u_coord','X','R_left','R_right','loops');
    fig = figure;
    set(fig,'visible','off');
    plot(X(:),u(:),'.-');
    titlename = sprintf('m = %d, J = %d, order = %d, lim = %.3d, N = %d', m, J, basis_order, type_limiter,ceil(loops/floor(loops/num_fig))+1);
    title(titlename);
    axis([R_left,R_right,-0.1,1.2]);
    print(fig,'-dpng',figname);
    err = -1;
end
if type_limiter > 0
    filename = sprintf('PME-m%.1f-c%.1f-p%.1f-basis%d-lim%.3d-J%d-ALL',m,c,p,basis_order,type_limiter,J);
    filename = sprintf('%s%s.mat',path_data,filename);
    save(filename,'m','c','p','R_left','R_right','dT', ...
        'J','basis_order','type_problem','type_limiter','mu', ...
        'loops','track_mean','track_osc','track_pos');
end
end

