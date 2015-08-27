clear;close all;

date = datestr(now,30);
path_data = sprintf('./data/ex_%s/',date);
path_report = sprintf('%sreport_%s.txt',path_data, date);
mkdir(path_data);
fid = fopen(path_report,'w');
fprintf(fid, 'Executed File: main_ex.m\nExecuted Time: %s\n\nOutput:\n',datestr(now));
fclose(fid);

m = 3;
R = 6;
dT = 0.5;
J = 50;
min_basis_order = 3;
max_basis_order = 4;
type_problem = 1;
list_limiter = [101,102,111,112,121,122,131,132];
err_table = zeros(2 * (max_basis_order - min_basis_order) + 2, length(list_limiter));
tic;
for i = 0: max_basis_order - min_basis_order
    basis_order = min_basis_order + i;
for index_lim = 1:length(list_limiter)
    type_limiter = list_limiter(index_lim);
        PARA = [m,0,0;
                -R,R,dT;
                J,basis_order,0;
                type_problem,type_limiter,0;];
        uI = @(x) BarenblattSolution(x,1,m);
        err = PME(PARA, uI, path_data, path_report);
        err_table(2 * i + 1, index_lim) = err;
        fid = fopen(path_report,'at');
        fprintf(fid, '%s | m = %d, J = %d, basis_order = %d, limiter = %d, err = %.4e\n', datestr(now), m, J, basis_order, type_limiter, err);
        fclose(fid);
end
end
t_all = toc;
fid = fopen(path_report,'at');
fprintf(fid, '\nFinish Time: %s\nTime used: %.4f', datestr(now),t_all);
fclose(fid);