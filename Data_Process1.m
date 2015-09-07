clear;
close all;
path_data = '.\data\ex_20150907T173652\';

num_moments = 200;
m = 3;
R = 6;
dT = 0.5;
J = 50;
min_basis_order = 1;
max_basis_order = 4;
type_problem = 1;
% list_limiter = [101,102,111,112,121,122,131,132];
list_limiter = [000,100,200];
for i = 0: max_basis_order - min_basis_order
    basis_order = min_basis_order + i;
for index_lim = 1:length(list_limiter)
    type_limiter = list_limiter(index_lim);
    filename = sprintf('%sPME-m%.1f-c%.1f-p%.1f-basis%d-lim%.3d-J%d-ALL.mat',path_data,m,0,0,basis_order,type_limiter,J);
    load(filename);
    figure;
    list_loop = 1:floor(loops/num_moments):loops;
    [x_mean,y_mean,value_mean] = find(track_mean(:,list_loop));
    [x_osc,y_osc,~] = find(track_osc(:,list_loop));
    [x_pos,y_pos,~] = find(track_pos(:,list_loop));
    k = 1/(m+1);
    interface_x = sqrt(2*m/k/(m-1)) * (list_loop / loops * dT + 1) .^ k;
    interface_x1 = ceil((interface_x + R) / (2*R/J));
    interface_x2 = ceil((-interface_x + R) / (2*R/J));
    interface_y = 1:length(list_loop);
    %plot(x_mean,y_mean,'s',x_osc,y_osc,'.',x_pos,y_pos,'o',interface_x1,interface_y,'r',interface_x2,interface_y,'r');
    plot(x_mean,y_mean,'s',interface_x1,interface_y,'r',interface_x2,interface_y,'r');
    axis([0,J+1,0,length(list_loop)+1]);
    if ~isempty(value_mean)
        value_mean = abs(value_mean);
        % fprintf('min(negative average) = %.4e, mean(negative average) = %.4e\n', max(value_mean), mean(value_mean));
        fprintf('min(negative average) = %.2e, time = %.2f\n', max(value_mean), time_consume);
    end
end
end