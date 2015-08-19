clear;close all;

date = datestr(now,30);
path_data = sprintf('./data/test_%s/',date);
path_report = sprintf('%sreport_%s.txt',path_data, date);
mkdir(path_data);
fid = fopen(path_report,'w');
fprintf(fid, 'Executed File: main_test.m\nExecuted Time: %s\n\nOutput:\n',datestr(now));
fclose(fid);

list_m = [3,5,7,9];
R = 6;
dT = 0.5;
J = 50;
min_basis_order = 1;
max_basis_order = 4;
type_problem = 0;
type_limiter = 0;
CFL_table = zeros(max_basis_order - min_basis_order + 1, length(list_m));
tic;
for index_m = 1:length(list_m)
    m = list_m(index_m);
    % The initial CFL will be enlarged if the assert in L_pme.m is wrong
    CFL = 10;
    CFL_step = 20;
    for i = 0: max_basis_order - min_basis_order
        basis_order = min_basis_order + i;
        err = 0/0;
        while isnan(err)
            try
                PARA = [m,0,0;
                        -R,R,dT;
                        J,basis_order,0;
                        type_problem,type_limiter,CFL;];
                uI = @(x) BarenblattSolution(x,1,m);
                err = PME(PARA, uI, path_data, path_report);
            catch ME
                CFL = CFL + CFL_step;
            end
        end
        CFL_table(i + 1, index_m) = CFL;
        fid = fopen(path_report,'at');
        fprintf(fid, '%s | m = %d, J = %d, basis_order = %d, limiter = %d\n', datestr(now), m, J, basis_order, type_limiter);
        fclose(fid);
    end
end
t_all = toc;
filename = sprintf('%s%s',path_data, 'CFL_table.mat');
save(filename,'CFL_table');
fid = fopen(path_report,'at');
fprintf(fid, '\nFinish Time: %s\nTime used: %.4f', datestr(now),t_all);
fclose(fid);