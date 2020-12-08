function plot_data3D_test(struct, order, slic, plt_color_cmds)
data = struct.data;
axes = struct.axes;
data_cell = fn_struct2cell(data);
axes_cell = fn_struct2cell(axes);

X = axes_cell{2,strcmpi(axes_cell(4,:), order{1})};
Y = axes_cell{2,strcmpi(axes_cell(4,:), order{2})};
Z = axes_cell{2,strcmpi(axes_cell(4,:), order{3})};
units = {'frequency','log_frequency' 'amplitude', 'offset';' [Hz]', ' [log(Hz)]', ' [V]', ' [V]'};

Xlbl = [strrep(order{1},'_','\_') units{2,strcmpi(units(1,:), order{1})}];
Ylbl = [strrep(order{2},'_','\_') units{2,strcmpi(units(1,:), order{2})}];
Zlbl = [strrep(order{3},'_','\_') units{2,strcmpi(units(1,:), order{3})}];

Xsl = slic{2,strcmpi(slic(1,:), order{1})};
Ysl = slic{2,strcmpi(slic(1,:), order{2})};
Zsl = slic{2,strcmpi(slic(1,:), order{3})};

figure
clf
[subrow, subcol] = subplot_min_rectangle(size(data_cell,2));
for i = 1:size(data_cell,2)
    subplot(subrow, subcol, i)
    slice(Y,X,Z, data_cell{2,i}, Ysl, Xsl, Zsl);
    grid on
    xlabel(Ylbl)
    ylabel(Xlbl)
    zlabel(Zlbl)
    title(data_cell{4,i})
    for c = plt_color_cmds
        eval(c{:});
    end
end
end

