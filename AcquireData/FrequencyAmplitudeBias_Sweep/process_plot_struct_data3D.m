function [fig, sbp, OutStruct, data_cell, axes_cell, NDims]  = process_plot_struct_data3D(struct, order, slic, varargin)
CharOrString = @(s) ischar(s) || isstring(s);
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('ax_range', {}, @iscell);
p.addParameter('val_range', {}, @iscell);
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
%% Limit by ax_range
xr = [-inf inf];
yr = [-inf inf];
zr = [-inf inf];    
if ~isempty(p.Results.ax_range)
    % Get upper and lower ax_range
    if any(strcmpi(p.Results.ax_range(1,:), order{1}))
        xr = p.Results.ax_range{2,strcmpi(p.Results.ax_range(1,:), order{1})};
    end
    if any(strcmpi(p.Results.ax_range(1,:), order{2}))
        yr = p.Results.ax_range{2,strcmpi(p.Results.ax_range(1,:), order{2})};
    end
    if any(strcmpi(p.Results.ax_range(1,:), order{3}))
        zr = p.Results.ax_range{2,strcmpi(p.Results.ax_range(1,:), order{3})}; 
    end
end
% Compute length in each dimension
xl = sum(X(:,1,1)>=xr(1) & X(:,1,1)<=xr(2));
yl = sum(Y(1,:,1)>=yr(1) & Y(1,:,1)<=yr(2));
zl = sum(Z(1,1,:)>=zr(1) & Z(1,1,:)<=zr(2));
% Pick elements within ax_range (inner box)
Xr = X>=xr(1) & X<=xr(2); 
Yr = Y>=yr(1) & Y<=yr(2);
Zr = Z>=zr(1) & Z<=zr(2);
Ir = Xr & Yr & Zr;
% reshape to 3D array
X = reshape(X(Ir), xl,yl,zl);
Y = reshape(Y(Ir), xl,yl,zl); 
Z = reshape(Z(Ir), xl,yl,zl);

%% Add axes vectors to struct
OutStruct.axesvec.(order{1}).vec = shiftdim(X(:,1,1));
OutStruct.axesvec.(order{2}).vec = shiftdim(Y(1,:,1));
OutStruct.axesvec.(order{3}).vec = shiftdim(Z(1,1,:));
% Addlabels
OutStruct.axesvec.(order{1}).label = Xlbl;
OutStruct.axesvec.(order{2}).label = Ylbl;
OutStruct.axesvec.(order{3}).label = Zlbl;


%% Slice planes
Xs = slic{2,strcmpi(slic(1,:), order{1})};
Xsl = Xs(Xs>=min(X,[],'all') & Xs<=max(X,[],'all'));
Ys = slic{2,strcmpi(slic(1,:), order{2})};
Ysl = Ys(Ys>=min(Y,[],'all') & Ys<=max(Y,[],'all'));
Zs = slic{2,strcmpi(slic(1,:), order{3})};
Zsl = Zs(Zs>=min(Z,[],'all') & Zs<=max(Z,[],'all'));


%% Process and limit by val_range
data_cell(1,:) = data_cell(4,:);
if ~isempty(p.Results.plt_select_data)
    data_cell = data_cell(:,cellfun(@(c) contains(c,p.Results.plt_select_data(1,:)),data_cell(1,:)));
end
for i = 1:size(data_cell,2)
    data_cell{2,i} = reshape(data_cell{2,i}(Ir), xl,yl,zl);
    if any(strcmpi(p.Results.plt_select_data(1,:),data_cell{1,i}))
        data_cell(3,i) = p.Results.plt_select_data(2,strcmpi(p.Results.plt_select_data(1,:),data_cell{1,i}));
        if any(strcmpi(p.Results.val_range(1,:),data_cell{1,i}))
            val_range = p.Results.val_range{2,strcmpi(p.Results.val_range(1,:),data_cell{1,i})};
            data_cell{2,i}(data_cell{2,i}<val_range(1) | data_cell{2,i}>val_range(2)) = NaN;
        end
    else
        data_cell(3,i) = data_cell(1,:);
    end
end
data_cell(4,:) = cellfun(@(c) c(1:(sum(find(c==' ', 1, 'last'))+isempty(find(c==' ', 1, 'last'))*length(c))-1), data_cell(3,:), 'UniformOutput', false);

%% 2D
NDims = 3;
Vec123 = 1:3;
Dim1 = size(X)==1;
FindDim1 = find(Dim1,1);
if FindDim1
    NDims = 2;
    PermVec = [Vec123(~Dim1) Vec123(Dim1)];
    X = permute(X, PermVec);
    Y = permute(Y, PermVec);
    Z = permute(Z, PermVec);
    for i = 1:size(data_cell,2)
        data_cell{2,i} = permute(data_cell{2,i}, PermVec);
    end
end
%% Out axes
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

if FindDim1
    axes_cell(:,:) = axes_cell(:,PermVec);
end
%% Plot
XYZ = {X,Y,Z};
fig = figure;
[subrow, subcol] = subplot_min_rectangle(size(data_cell,2));
UserDataStruct.axesvec = OutStruct.axesvec;
for i = 1:size(data_cell,2)
    UserDataStruct.type = data_cell{3,i};
    UserDataStruct.data = data_cell{2,i};
    s =  subplot(subrow, subcol, i);
    if FindDim1
        surf(XYZ{~Dim1}, data_cell{2,i})
        xlabel(axes_cell{3,1})
        ylabel(axes_cell{3,2})
    else
        slice(XYZ{[2 1 3]}, data_cell{2,i}, Ysl, Xsl, Zsl);
        xlabel(Ylbl)
        ylabel(Xlbl)
        zlabel(Zlbl)
    end
    s.UserData = UserDataStruct;
    s.Title.UserData = data_cell{4,i};
    sbp(ceil(i/subcol), mod(i,subcol) + (mod(i,subcol)==0)*subcol) = s;
    if ~isempty(p.Results.title)
        title([p.Results.title ' ' data_cell{3,i}])
    else
        title(data_cell{3,i})
    end
    for c = p.Results.plt_cmds(2,strcmpi(p.Results.plt_cmds(1,:),data_cell{1,i}) | cellfun(@(c) isempty(c),p.Results.plt_cmds(1,:)))
        eval(c{:});
    end
end
OutStruct.data = data_cell;
OutStruct.axes = axes_cell;
set(fig,'Position',get(0,'Screensize'));
end

