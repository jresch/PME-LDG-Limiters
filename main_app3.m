clear;close all;
% main_app3.m
%   Objective:  Stores the moments of the numerical solution.
%   Output:     Screenshots of the solution.
%   Problem:    PME with absorption.
%   Initial:    |Sin(x)| with platform.
%   Solution:   Unknown.
% PARA = [  m, c, p;
%           R_left, R_right, dT;
%           J, basis_order, 0;
%           type_problem, type_limiter, CFL;];
date = datestr(now,30);
path_data = sprintf('./data/app3_%s/',date);
path_report = sprintf('%sreport_%s.txt',path_data, date);
mkdir(path_data);
fid = fopen(path_report,'w');
fprintf(fid, 'Executed File: main_app3.m\nExecuted Time: %s\n\nOutput:\n',datestr(now));
fclose(fid);

list_mcp = [1.5,1,0.1;
            1.9,1,0.1;
            2.5,1,0.1;
            2.5,1,0.5;
            2.5,0,0.1;
            2.5,10,0.1;];
R = 6;
dT = 1;
J = 50;
min_basis_order = 1;
max_basis_order = 1;
type_problem = 2;
type_limiter = 1;
tic;
for index_m = 1:size(list_mcp,1)
    m = list_mcp(index_m, 1);
    c = list_mcp(index_m, 2);
    p = list_mcp(index_m, 3);
    for i = 0: max_basis_order - min_basis_order
        basis_order = min_basis_order + i;
        PARA = [m,c,p;
                -R,R,dT;
                J,basis_order,0;
                type_problem,type_limiter,0;];
        uI = @(x) (abs(x)>=pi/6) .* (abs(x)<=pi) .* abs(sin(x)) + (abs(x)<pi/6) * 0.5;
        path_dir = sprintf('%s/m%.1f-c%.1f-p%.1f/',path_data,m,c,p);mkdir(path_dir);
        PME(PARA, uI, path_dir, path_report);
        fid = fopen(path_report,'at');
        fprintf(fid, '%s | m%.1f, c%.1f, p%.1f, J = %d, basis_order = %d\n', datestr(now), m, c, p, J, basis_order);
        fclose(fid);
    end
end
t_all = toc;
fid = fopen(path_report,'at');
fprintf(fid, '\nFinish Time: %s\nTime used: %.4f', datestr(now),t_all);
fclose(fid);