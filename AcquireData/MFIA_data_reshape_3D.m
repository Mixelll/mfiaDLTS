function [data_desired_order_3D, data_sweep_order_3D] = MFIA_data_reshape_3D(data_sweep_order_struct_vec, desired_order, sweep_order, sweep_pts, frequency_vec, amplitude_vec, offset_vec)
data_desired_order_3D = [];
data_sweep_order_3D = [];
pos = cellfun(@(c) find(strcmp(c, sweep_order)), desired_order);
data_in1 = data_sweep_order_struct_vec(1);
read_param_cell = fn_struct2cell(data_in1);

for c = read_param_cell(3,:)
    eval(['data_sweep_order_3D' c{:} '= reshape(cell2mat(cellfun(@transpose, {data_sweep_order_struct_vec' c{:} '})), sweep_pts, length(' sweep_order{2} '_vec), length(' sweep_order{3} '_vec));'])
    eval(['data_desired_order_3D' c{:} '= permute(data_desired_order_3D' c{:} ',' pos ');'])
end    
end

