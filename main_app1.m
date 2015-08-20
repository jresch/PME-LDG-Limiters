clear;close all;

date = datestr(now,30);
path_data = sprintf('./data/app1_%s/',date);
path_report = sprintf('%sreport_%s.txt',path_data, date);
mkdir(path_data);
fid = fopen(path_report,'w');
fprintf(fid, 'Executed File: main_app1.m\nExecuted Time: %s\n\nOutput:\n',datestr(now));
fclose(fid);

num_refine = 0;
list_m = [3];
R = 6;
dT = 3;
J = 50;
min_basis_order = 3;
max_basis_order = 3;
type_problem = 2;
type_limiter = 3;
tic;
for index_m = 1:length(list_m)
    m = list_m(index_m);
    for i = 0: max_basis_order - min_basis_order
        basis_order = min_basis_order + i;
        PARA = [m,0,0;
                -R,R,dT;
                J,basis_order,0;
                type_problem,type_limiter,0;];
        uI = @(x) (abs(x)>=1) .* (abs(x)<=3);
        path_data1 = sprintf('%s/Box1/',path_data);mkdir(path_data1);
        PME(PARA, uI, path_data1, path_report);
        fid = fopen(path_report,'at');
        fprintf(fid, '%s | Box1| m = %d, J = %d, basis_order = %d, limiter = %d\n', datestr(now), m, J, basis_order, type_limiter);
        fclose(fid);
        uI = @(x) 0.5 * ((abs(x)>=1) .* (abs(x)<=3) + (x>=1) .* (x<=3));
        path_data2 = sprintf('%s/Box2/',path_data);mkdir(path_data2);
        PME(PARA, uI, path_data2, path_report);
        fid = fopen(path_report,'at');
        fprintf(fid, '%s | Box2| m = %d, J = %d, basis_order = %d, limiter = %d\n', datestr(now), m, J, basis_order, type_limiter);
        fclose(fid);
    end
end
t_all = toc;
fid = fopen(path_report,'at');
fprintf(fid, '\nFinish Time: %s\nTime used: %.4f', datestr(now),t_all);
fclose(fid);