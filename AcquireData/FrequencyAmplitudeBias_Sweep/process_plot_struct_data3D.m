function [fig, sbp, data_cell, axes_cell]  = process_plot_struct_data3D(struct, order, slic, varargin)
CharOrString = @(s) ischar(s) || isstring(s);
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('range', {}, @iscell);
p.addParameter('omit', {}, @iscell);
p.addParameter('plt_select_data', {}, @iscell);
p.addParameter('plt_cmds', {}, @iscell);
p.addParameter('title', '', CharOrString);
p.parse(varargin{:});

data = struct.data;
axes = struct.axes;
data_cell = fn_struct2cell(data);
axes_cell = fn_struct2cell(axes);
% Frequency, amplitude and offset grid
AxO1 = strcmpi(axes_cell(4,:), order{1});
AxO2 = strcmpi(axes_cell(4,:), order{2});
AxO3 = strcmpi(axes_cell(4,:), order{3});
X = axes_cell{2,AxO1};
Y = axes_cell{2,AxO2};
Z = axes_cell{2,AxO3};
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
if ~isempty(p.Results.range)
    % Get upper and lower range
    xr = p.Results.range{2,strcmpi(p.Results.range(1,:), order{1})};
    yr = p.Results.range{2,strcmpi(p.Results.range(1,:), order{2})};
    zr = p.Results.range{2,strcmpi(p.Results.range(1,:), order{3})};
else
    xr = [-inf inf];
    yr = [-inf inf];
    zr = [-inf inf];
end
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

axes_cell{1,1} = order{1};
axes_cell{1,2} = order{2};
axes_cell{1,3} = order{3};
axes_cell{2,1} = X;
axes_cell{2,2} = Y;
axes_cell{2,3} = Z;
axes_cell{3,1} = Xlbl;
axes_cell{3,2} = Ylbl;
axes_cell{3,3} = Zlbl;
axes_cell = axes_cell(1:3,1:3);
axes_cell(4,:) = axes_cell(1,:);
%% Process
data_cell(1,:) = data_cell(4,:);
for i = 1:size(data_cell,2)
    data_cell{2,i} = reshape(data_cell{2,i}(Ir), xl,yl,zl);
    if ~isempty(p.Results.plt_select_data(2,strcmpi(p.Results.plt_select_data(1,:),data_cell{1,i})))
        data_cell(3,i) = p.Results.plt_select_data(2,strcmpi(p.Results.plt_select_data(1,:),data_cell{1,i}));
    else
        data_cell(3,i) = data_cell(1,:);
    end
end
data_cell(4,:) = cellfun(@(c) c(1:(sum(find(c==' ', 1, 'last'))+isempty(find(c==' ', 1, 'last'))*length(c))-1), data_cell(3,:), 'UniformOutput', false);

%% Plot
fig = figure;
if ~isempty(p.Results.plt_select_data)
    data_cell = data_cell(:,cellfun(@(c) contains(c,p.Results.plt_select_data(1,:)),data_cell(1,:)));
end
[subrow, subcol] = subplot_min_rectangle(size(data_cell,2));

for i = 1:size(data_cell,2)
    s =  subplot(subrow, subcol, i);
    sbp(subrow, subcol) = s;
    slice(Y,X,Z, data_cell{2,i}, Ysl, Xsl, Zsl);
    xlabel(Ylbl)
    ylabel(Xlbl)
    zlabel(Zlbl)
    if ~isempty(p.Results.title)
        title([p.Results.title ' ' data_cell{3,i}])
    else
        title(data_cell{3,i})
    end
    for c = p.Results.plt_cmds(2,strcmpi(p.Results.plt_cmds(1,:),data_cell{1,i}) | cellfun(@(c) isempty(c),p.Results.plt_cmds(1,:)))
        eval(c{:});
    end
end
set(fig,'Position',get(0,'Screensize'));
end

