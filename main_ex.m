clear;close all;

date = datestr(now,30);
path_data = sprintf('./data/ex_%s/',date);
path_report = sprintf('%sreport_%s.txt',path_data, date);
mkdir(path_data);
fid = fopen(path_report,'w');
fprintf(fid, 'Executed File: main_ex.m\nExecuted Time: %s\n\nOutput:\n',datestr(now));
fclose(fid);

num_refine = 0;
list_m = [4,5];
R = 6;
dT = 0.5;
J_initial = 50;
min_basis_order = 1;
max_basis_order = 4;
type_problem = 1;
type_limiter = 3;
err_table = zeros(num_refine + 1, 2 * (max_basis_order - min_basis_order) + 2, length(list_m));
tic;
for index_m = 1:length(list_m)
    m = list_m(index_m);
for i = 0: max_basis_order - min_basis_order
    basis_order = min_basis_order + i;
    for j = 0:num_refine
        PARA = [m,0,0;
                -R,R,dT;
                J_initial*2^j,basis_order,0;
                type_problem,type_limiter,0;];
            uI = @(x) BarenblattSolution(x,1,m);
        err_table(j + 1, 2 * i + 1, index_m) = PME(PARA, uI, path_data, path_report);
        fid = fopen(path_report,'at');
        fprintf(fid, '%s | m = %d, J = %d, basis_order = %d, limiter = %d\n', datestr(now), m, J_initial * 2 ^ j, basis_order, type_limiter);
        fclose(fid);
    end
    err_table(1, 2 * i + 2, index_m) = 0;
    err_table(2:end, 2 * i + 2, index_m) = log2(err_table(1:end-1,2*i+1,index_m)./err_table(2:end,2*i+1,index_m));
end
end
t_all = toc;
fid = fopen(path_report,'at');
fprintf(fid, '\nFinish Time: %s\nTime used: %.4f', datestr(now),t_all);
fclose(fid);