path_data = '.\data\app1_20150819T144453\Box1\';

num_moments = 25;
list_m = [3];
J = 50;
min_basis_order = 2;
max_basis_order = 2;
type_limiter = 1;
for index_m = 1:length(list_m)
    m = list_m(index_m);
    for i = 0: max_basis_order - min_basis_order
        basis_order = min_basis_order + i;
        filename = sprintf('%sPME-m%.1f-c%.1f-p%.1f-basis%d-lim%d-J%d-ALL.mat',path_data,m,0,0,basis_order,type_limiter,J);
        load(filename);
        figure;
        list_loop = 1:floor(loops/num_moments):loops;
        [x_mean,y_mean,value_mean] = find(track_mean(:,list_loop));
        [x_osc,y_osc,~] = find(track_osc(:,list_loop));
        [x_pos,y_pos,~] = find(track_pos(:,list_loop));
        plot(x_mean,y_mean,'s',x_osc,y_osc,'.',x_pos,y_pos,'o');
    end
end