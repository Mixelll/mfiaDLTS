function fn_plot_data3D(struct, order, range, slic, plt_select_data, plt_cmds)
data = struct.data;
axes = struct.axes;
data_cell = fn_struct2cell(data);
axes_cell = fn_struct2cell(axes);
% Frequency, amplitude and offset grid 
X = axes_cell{2,strcmpi(axes_cell(4,:), order{1})};
Y = axes_cell{2,strcmpi(axes_cell(4,:), order{2})};
Z = axes_cell{2,strcmpi(axes_cell(4,:), order{3})};
% Axis labels
axes_lbls = {'frequency','log_frequency' 'amplitude', 'offset';'Test Signal Frequency [Hz]', 'Test Signal Frequency [log(Hz)]', 'Test Signal Amplitude [V]', 'Bias on Si [V]'};
Xlbl = axes_lbls{2,strcmpi(axes_lbls(1,:), order{1})};
Ylbl = axes_lbls{2,strcmpi(axes_lbls(1,:), order{2})};
Zlbl = axes_lbls{2,strcmpi(axes_lbls(1,:), order{3})};
% Slice planes
Xsl = slic{2,strcmpi(slic(1,:), order{1})};
Ysl = slic{2,strcmpi(slic(1,:), order{2})};
Zsl = slic{2,strcmpi(slic(1,:), order{3})};
%% Plot range
% Get upper and lower range
xr = range{2,strcmpi(range(1,:), order{1})};
yr = range{2,strcmpi(range(1,:), order{2})};
zr = range{2,strcmpi(range(1,:), order{3})};
% Compute length in each dimension
xl = sum(X(:,1,1)>=xr(1) & X(:,1,1)<=xr(2));
yl = sum(Y(1,:,1)>=yr(1) & Y(1,:,1)<=yr(2));
zl = sum(Z(1,1,:)>=zr(1) & Z(1,1,:)<=zr(2));
% Pick elements within range (inner box)
Xr = X>=xr(1) & X<=xr(2); 
Yr = Y>=yr(1) & Y<=yr(2);
Zr = Z>=zr(1) & Z<=zr(2);
Ir = Xr & Yr & Zr;
% reshape to 3D array
X = reshape(X(Ir), xl,yl,zl);
Y = reshape(Y(Ir), xl,yl,zl); 
Z = reshape(Z(Ir), xl,yl,zl); 
%%
% plot
figure
clf
if ~isempty(plt_select_data)
    data_cell = data_cell(:,cellfun(@(c) contains(c,plt_select_data(1,:)),data_cell(1,:)));
end
[subrow, subcol] = subplot_min_rectangle(size(data_cell,2));
for i = 1:size(data_cell,2)
    subplot(subrow, subcol, i)
    slice(Y,X,Z, reshape(data_cell{2,i}(Ir), xl,yl,zl), Ysl, Xsl, Zsl);
    xlabel(Ylbl)
    ylabel(Xlbl)
    zlabel(Zlbl)
    title(plt_select_data(2,strcmpi(plt_select_data(1,:),data_cell{4,i})))
    for c = plt_cmds
        eval(c{:});
    end
end
end

