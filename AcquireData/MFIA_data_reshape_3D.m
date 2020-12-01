function [data_desired_order_3D, data_sweep_order_3D] = MFIA_data_reshape_3D(data_sweep_order_struct_vec, desired_order, sweep_order, pts_frequency, pts_amplitude, pts_offset, frequency_vec, amplitude_vec, offset_vec)
data_desired_order_3D = [];
data_sweep_order_3D = [];
pos = cellfun(@(c) find(strcmpi(c, sweep_order)), desired_order);
data_in1 = data_sweep_order_struct_vec(1);
read_param_cell = fn_struct2cell(data_in1);
read_param_cell = read_param_cell(:,~contains(read_param_cell(1,:), sweep_order(2:3),'IgnoreCase',true));
sweep_pts = eval(['pts_' sweep_order{1}]);
l2 = eval(['pts_' sweep_order{2}]);
l3 = eval(['pts_' sweep_order{3}]);
for c = read_param_cell(3,:)
    tmp = nan(sweep_pts,l2*l3);
    tmpperm = [];
    for i = 1:size(tmp,2)
        tmp(:,i) = eval(['data_sweep_order_struct_vec(i)' c{:} '.''']);
    end
    eval(['data_sweep_order_3D' c{:} '= reshape(tmp, sweep_pts, l2, l3);'])
    tmpperm = permute(eval(['data_sweep_order_3D' c{:}]), pos);
    eval(['data_desired_order_3D' c{:} '= tmpperm;'])
end    
end

