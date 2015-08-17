function [ err ] = PME( PARA, uI, path_data, path_report )
% PME.m
%   This function is the solver of PME and it returns the error.
%   When type_problem = 0 or 1, it returns the error,
%   otherwise it returns -1.
%   There are three global variables which store the information of
%   limiters.
%       count_mean  where the mean of u is < 0.
%       count_osc   where the oscillation exist.
%       count_pos   where the negative value exist.

global m;
global c;
global p;
global weights;
global dx;
global psi;
global psi_z;
global mu;
global count_mean;
global count_osc;
global count_pos;

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
% Gauss-Lobatto Quadrature
[points,weights] = Quadrature_Set();
% Mesh
[grid,X] = Mesh_Set(points,R_left,R_right,J);
% Time Step
dx = max(grid(2,:)) * 2;
if type_problem == 0
    CFL =  PARA(4,3);
    loops = ceil(dT / dx / dx * CFL);
else
    CFL = [ 50, 70, 100, 130;
            170, 270, 390, 490;
            490, 790, 1110, 1450;
            1110, 1850, 2610, 3370;];
    loops = ceil(dT / dx / dx * CFL(basis_order,floor(m/2)));
end
dt = dT / loops;
% Basis
[psi,psi_z] = Basis_Set(points,basis_order);
% Limiters
count_mean = sparse(J,loops);
count_osc = sparse(J,loops);
count_pos = sparse(J,loops);
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
    u_coord = L_pme(u0_coord, l);
    u1_coord = u0_coord + u_coord * dt;
    u1_coord = Limiter(u1_coord,l,type_limiter);
    % step 2
    u_coord = L_pme(u1_coord, l);
    u2_coord = u0_coord * 3/4 + u1_coord /4 + u_coord * dt/4;
    u2_coord = Limiter(u2_coord,l,type_limiter);
    % step 3
    u_coord = L_pme(u2_coord, l);
    u3_coord = u0_coord / 3 + u2_coord * 2/3 + u_coord * dt*2/3;
    u3_coord = Limiter(u3_coord,l,type_limiter);
    u_coord = u3_coord;
    if (type_problem == 2) && (mod(l,floor(loops/100)) == 1)
        u = psi * u_coord;
        filename = sprintf('PME-m%.1f-c%.1f-p%.1f-basis%d-lim%d-J%d-N%d',m,c,p,basis_order,type_limiter,J,ceil(l/floor(loops/100)));
        figname = sprintf('%s%s.png',path_data,filename);
        filename = sprintf('%s%s.mat',path_data,filename);
        save(filename,'u','u_coord','X','R_left','R_right','loops');
        fig = figure;
        set(fig,'visible','off');
        plot(X(:),u(:),'.-');
        axis([R_left,R_right,-0.1,1.2]);
        print(fig,'-dpng',figname);
    end
end
u = psi * u_coord;
if type_problem == 0
    u_exact_pme = BarenblattSolution(X, 1+dT, m);
    err = sqrt(dx/2 * sum(weights * (u_exact_pme - psi * u_coord) .^ 2));
    titlename = sprintf('m = %d, J = %d, order = %d', m, J, basis_order);
    figure;plot(X(:),u(:),'b',X(:),u_exact_pme(:),'r');title(titlename);
elseif type_problem == 1
    u_exact_pme = BarenblattSolution(X, 1+dT, m);
    err = sqrt(dx/2 * sum(weights * (u_exact_pme - psi * u_coord) .^ 2));
    titlename = sprintf('m = %d, J = %d, order = %d', m, J, basis_order);
    figure;plot(X(:),u(:),'b',X(:),u_exact_pme(:),'r');title(titlename);
elseif type_problem == 2
    filename = sprintf('PME-m%.1f-c%.1f-p%.1f-basis%d-lim%d-J%d-N%d',m,c,p,basis_order,type_limiter,J,ceil(loops/floor(loops/100))+1);
    figname = sprintf('%s%s.png',path_data,filename);
    filename = sprintf('%s%s.mat',path_data,filename);
    save(filename,'u','u_coord','X','R_left','R_right','loops');
    fig = figure;
    set(fig,'visible','off');
    plot(X(:),u(:),'.-');
    axis([R_left,R_right,-0.1,1.2]);
    print(fig,'-dpng',figname);
    err = -1;
%     interface_func = @(t) sqrt(2*m*(m+1)/(m-1)) * t.^(1/(m+1));
%     interface_cell = [  ceil((-interface_func((1:loops)*dt+1)-R_left)/dx),...
%                         ceil((interface_func((1:loops)*dt+1)-R_left)/dx)];
%     interface = sparse(interface_cell,[1:loops,1:loops],[1:loops,1:loops],J,loops);
%     figure;
%     subplot(1,2,1);plot(X(:),u(:),'b',X(:),u_exact_pme(:),'r');title(titlename);
%     subplot(1,2,2);
%     step = 15;
%     [x0,y0,~] = find(interface(:,1:floor(end/step):end));
%     [x1,y1,z1] = find(count_pos_diff(:,1:floor(end/step):end));
%     [x2,y2,z2] = find(count_pos_side(:,1:floor(end/step):end));
%     [x3,y3,z3] = find(count_pos_mean(:,1:floor(end/step):end));
%     plot(x0,y0,'s',x1,y1,'o',x2,y2,'+',x3,y3,'*');
%     legend('true','diff','side','mean');
%     xlim([0,J+1]);set(gca,'xtick',[]);set(gca,'ytick',[]);
end

end
