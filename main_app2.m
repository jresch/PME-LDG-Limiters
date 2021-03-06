clear;close all;

date = datestr(now,30);
path_data = sprintf('./data/app2_%s/',date);
path_report = sprintf('%sreport_%s.txt',path_data, date);
mkdir(path_data);
fid = fopen(path_report,'w');
fprintf(fid, 'Executed File: main_app2.m\nExecuted Time: %s\n\nOutput:\n',datestr(now));
fclose(fid);

list_m = [3];
R = 3;
dT = 3;
J = 25;
min_basis_order = 1;
max_basis_order = 1;
type_problem = 2;
type_limiter = 000;
tic;
for index_m = 1:length(list_m)
    m = list_m(index_m);
    for i = 0: max_basis_order - min_basis_order
        basis_order = min_basis_order + i;
        PARA = [m,0,0;
                -R,R,dT;
                J,basis_order,0;
                type_problem,type_limiter,0;];
        uI = @(x) cos(x) .* (abs(x)<=pi/2);
        PME(PARA, uI, path_data, path_report);
        fid = fopen(path_report,'at');
        fprintf(fid, '%s | m = %d, J = %d, basis_order = %d, limiter = %d\n', datestr(now), m, J, basis_order, type_limiter);
        fclose(fid);
    end
end
t_all = toc;
fid = fopen(path_report,'at');
fprintf(fid, '\nFinish Time: %s\nTime used: %.4f', datestr(now),t_all);
fclose(fid);