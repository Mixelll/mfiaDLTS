function [out_desired_order_3D, out_sweep_order_3D] = MFIA_data_reshape_3D(data_sweep_order_struct_vec, desired_order, sweep_order, pts_frequency, pts_amplitude, pts_offset, frequency_vec, amplitude_vec, offset_vec)
data_desired_order_3D = [];
data_sweep_order_3D = [];
axes_desired_order_3D = [];
axes_sweep_order_3D = [];
acaxes_desired_order_3D = [];
acaxes_sweep_order_3D = [];

out_desired_order_3D = [];
out_sweep_order_3D = [];

pos = cellfun(@(c) find(strcmpi(c, sweep_order)),desired_order);
data_in1 = data_sweep_order_struct_vec(1);
read_param_cell = fn_struct2cell(data_in1);
read_param_cell_val = read_param_cell(:,~contains(read_param_cell(1,:),sweep_order(2:3),'IgnoreCase',true));
read_param_cell_ax = read_param_cell(:,contains(read_param_cell(1,:),sweep_order(2:3),'IgnoreCase',true) & ~contains(read_param_cell(1,:),'actual','IgnoreCase',true));
read_param_cell_acax = read_param_cell(:,contains(read_param_cell(1,:),sweep_order(2:3),'IgnoreCase',true) & contains(read_param_cell(1,:),'actual','IgnoreCase',true));
sweep_pts = eval(['pts_' sweep_order{1}]);
l2 = eval(['pts_' sweep_order{2}]);
l3 = eval(['pts_' sweep_order{3}]);

%set values (axes)
for c = read_param_cell_ax
    tmp = nan(sweep_pts,l2*l3);
    for i = 1:size(tmp,2)
        tmp(:,i) = repmat(eval(['data_sweep_order_struct_vec(i)' c{3}]), sweep_pts, 1);
    end
    
    eval(['axes_sweep_order_3D.' c{4} '= reshape(tmp, sweep_pts, l2, l3);'])
    tmpperm = permute(eval(['axes_sweep_order_3D.' c{4}]), pos);
    eval(['axes_desired_order_3D.' c{4} '= tmpperm;'])
end

%actual values (acaxes)
for c = read_param_cell_acax
    tmp = nan(sweep_pts,l2*l3);
    for i = 1:size(tmp,2)
        tmp(:,i) = repmat(eval(['data_sweep_order_struct_vec(i)' c{3}]), sweep_pts, 1);
    end
    
    eval(['acaxes_sweep_order_3D.' c{4} '= reshape(tmp, sweep_pts, l2, l3);'])
    tmpperm = permute(eval(['acaxes_sweep_order_3D.' c{4}]), pos);
    eval(['acaxes_desired_order_3D.' c{4} '= tmpperm;'])
end

% data (measured values)
for c = read_param_cell_val(3,:)
    tmp = nan(sweep_pts,l2*l3);
    for i = 1:size(tmp,2)
        tmp(:,i) = eval(['data_sweep_order_struct_vec(i)' c{:} '.''']);
    end
    
    eval(['data_sweep_order_3D' c{:} '= reshape(tmp, sweep_pts, l2, l3);'])
    tmpperm = permute(eval(['data_sweep_order_3D' c{:}]), pos);
    eval(['data_desired_order_3D' c{:} '= tmpperm;'])
end
    
% sweeper grid
cn = read_param_cell(:,contains(read_param_cell(1,:),'grid','IgnoreCase',true));
c = cn(:,1);

eval(['axes_desired_order_3D.' sweep_order{1} '=data_desired_order_3D' c{3,:} ';'])
eval(['axes_sweep_order_3D.' sweep_order{1} '=data_sweep_order_3D' c{3,:} ';'])
eval(['acaxes_desired_order_3D.' sweep_order{1} '=axes_desired_order_3D.' sweep_order{1} ';'])
eval(['acaxes_sweep_order_3D.' sweep_order{1} '=axes_sweep_order_3D.' sweep_order{1} ';'])

% removes grid (X axis) from measured data
for c = cn
    grid_field = strrep(c{3},'.grid','');
    v = rmfield(eval(['data_desired_order_3D' grid_field]), 'grid'); 
    k = rmfield(eval(['data_sweep_order_3D' grid_field]), 'grid');
    if isempty(fieldnames(v))
        data_desired_order_3D = rmfield(data_desired_order_3D, grid_field(2:end));
        data_sweep_order_3D = rmfield(data_sweep_order_3D, grid_field(2:end));
    else
        eval(['data_desired_order_3D' grid_field '=v;'])
        eval(['data_sweep_order_3D' grid_field '=k;'])
    end
end

% give frequency in log10 values
axes_desired_order_3D.log_frequency = log10(axes_desired_order_3D.frequency);
axes_sweep_order_3D.log_frequency = log10(axes_sweep_order_3D.frequency);
acaxes_desired_order_3D.log_frequency = log10(acaxes_desired_order_3D.frequency);
acaxes_sweep_order_3D.log_frequency = log10(acaxes_sweep_order_3D.frequency);

out_desired_order_3D.data = data_desired_order_3D;
out_desired_order_3D.axes = axes_desired_order_3D;
out_desired_order_3D.acaxes = acaxes_desired_order_3D;
out_sweep_order_3D.data = data_sweep_order_3D;
out_sweep_order_3D.axes = axes_sweep_order_3D;
out_sweep_order_3D.acaxes = acaxes_sweep_order_3D;
end

