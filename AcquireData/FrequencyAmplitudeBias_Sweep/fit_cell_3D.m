function [fig, sbp, OutStructArray, OutAppendPath, StrFitUscore] = fit_cell_3D(DataCell, AxesCell, fit_func, fit_axis, varargin)
CharOrString = @(s) ischar(s) || isstring(s);
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('val_range', {}, @iscell);
p.addParameter('plt_select_data', {}, @iscell);
p.addParameter('plt_cmds', {}, @iscell);
p.addParameter('title', '', CharOrString);
p.addParameter('plot_fit', []);
p.addParameter('str', '', CharOrString);
p.parse(varargin{:});

if ~isempty(p.Results.plot_fit)
    if isfield(p.Results.plot_fit, 'func')
        plot_func = p.Results.plot_fit.func;
    else
        plot_func = @plot;
    end
    if isfield(p.Results.plot_fit, 'visible')
        plot_visible = p.Results.plot_fit.visible;
    else
        plot_visible = false;
    end
    if isfield(p.Results.plot_fit, 'properties')
        plot_properties = p.Results.plot_fit.properties;
    else
        plot_properties = {};
    end
    if isfield(p.Results.plot_fit, 'commands')
        plot_commands = p.Results.plot_fit.commands;
    else
        plot_commands = {};
    end
    if isfield(p.Results.plot_fit, 'screensize')
        plot_screensize = p.Results.plot_fit.screensize;
    else
        plot_screensize = 0;
    end
end

StartTime = datestr(now, 'yyyy-mm-dd HH-MM-SS');

d123 = 1:3;
ftAxPos = strcmpi(AxesCell(1,:),fit_axis);
perm_vec = [d123(ftAxPos) d123(~ftAxPos)];
AxesCell(2,:) = cellfun(@(c) permute(c,perm_vec), AxesCell(2,:) , 'UniformOutput',false);
AxesCell(:,:) = AxesCell(:,perm_vec);
DataCell(2,:) = cellfun(@(c) permute(c,perm_vec), DataCell(2,:) , 'UniformOutput',false);

% Create axes matrices
FitAxMat = AxesCell{2,1};
AxJ = shiftdim(mean(AxesCell{2,2}));
AxK = shiftdim(mean(AxesCell{2,3}));

% create axes vectors struct
axesvec.(AxesCell{1,2}).vec = shiftdim(AxJ(:,1,1));
axesvec.(AxesCell{1,3}).vec = shiftdim(AxK(1,:,1));
% Addlabels
axesvec.(AxesCell{1,2}).label = AxesCell{3,2};
axesvec.(AxesCell{1,3}).label = AxesCell{3,3};

if ~isempty(p.Results.plt_select_data)
    DataCell = DataCell(:,cellfun(@(c) contains(c,p.Results.plt_select_data(1,:),'IgnoreCase',true), DataCell(1,:))...
        | cellfun(@(c) contains(c,p.Results.plt_select_data(1,:),'IgnoreCase',true), DataCell(3,:)));
end

cd_mat = DataCell{2,1};
[~, ~, ~, fit_cell0] = fit_func(FitAxMat(:,round(size(FitAxMat,2)/2),round(size(FitAxMat,3)/2)), cd_mat(:,round(size(FitAxMat,2)/2),round(size(FitAxMat,3)/2)));
fmax = size(fit_cell0,2);
FitOut = zeros([size(FitAxMat, [2 3]) fmax]);
StrFitUscore = [];
for c = fit_cell0(1,:)
    StrFitUscore = [StrFitUscore '_' c{:}];
end
StrFitUscore(1) = '';
FitTitles = fit_cell0(3,:);
[subrow, subcol] = subplot_min_rectangle(size(DataCell,2)*fmax);
fig = figure;
i=1;
OutStructArray = [];
for cd = DataCell
    cd_mat = cd{2};
    for j=1:size(cd_mat, 2)
        for k=1:size(cd_mat, 3)
            [fit_str, OutFitFunc, span, fit_cell] = fit_func(FitAxMat(:,j,k), cd_mat(:,j,k));
            fit_progress = ['j-k=' num2str(j) '-' num2str(k) ' ' AxesCell{4,2} '=' num2str(AxJ(j,1),3) ' ' AxesCell{4,3} '=' num2str(AxK(j,k),3)];
            fit_progress_title = ['j-k=' num2str(j) '-' num2str(k) newline AxesCell{4,2} '=' num2str(AxJ(j,1),3) ' ' AxesCell{4,3} '=' num2str(AxK(j,k),3)];
            if contains(fit_str, 'complex','IgnoreCase',true)
                disp(fit_progress)
            end
            if ~isempty(p.Results.plot_fit)
                ax = s_plot(FitAxMat(:,j,k), cd_mat(:,j,k), plot_properties, plot_commands, [p.Results.title ' ' cd{4} newline strrep(fit_progress_title,'_','\_')], fit_str, '', AxesCell{3,1}, cd{3}, 15, plot_func, ~plot_visible, plot_screensize);
                if ~isempty(OutFitFunc)
                    hold(ax(1),'on');
                    fplot(OutFitFunc, span)
                    hold(ax(1),'off');
                end
                if isfield(p.Results.plot_fit, 'savepath')
                    SavePath = [p.Results.plot_fit.savepath '\' cd{4} '-' fit_axis '\' StrFitUscore ' ' p.Results.str StartTime];
                    if ~isempty(OutFitFunc)
                        RealSavePath = [SavePath '\success' ];
                    else
                        RealSavePath = [SavePath '\fail' ];
                    end
                    if ~exist(RealSavePath, 'dir')
                        mkdir(RealSavePath)
                    end
                    saveas(ax(1), [RealSavePath '\' fit_progress '.png']);
                end
            end
            for f=1:fmax
                FitOut(j,k,f) = fit_cell{2,f};
            end            
        end
    end
    OutStruct.name = [cd{4} '_' fit_axis '_' StrFitUscore];
    OutStruct.title = [cd{3} '-' fit_axis ' ' strrep(StrFitUscore,'_', ', ')];
    OutFitCellInner = cell(3,f);
    for f=1:fmax
        RangedFitOut = FitOut(:,:,f);
        RangedFitPlot = FitOut(:,:,f);
        if any(strcmpi(p.Results.val_range(1,:),fit_cell{1,f}))
        	val_range = p.Results.val_range{2,strcmpi(p.Results.val_range(1,:),fit_cell{1,f})};
            RangedFitPlot(RangedFitPlot<val_range(1) | RangedFitPlot>val_range(2)) = NaN;
%             eval(['OutStruct.metadata.' fit_cell{1,f} '_range = val_range;'])
        end
        figure(fig);
        UserDataStruct.axesvec = axesvec;
        UserDataStruct.type = FitTitles{f};
        s = subplot(subrow, subcol, i);
        surf(AxJ, AxK, RangedFitPlot);
        s.UserData = UserDataStruct;
        s.Title.UserData = FitTitles{f};
        sbp(ceil(i/subcol), mod(i,subcol) + (mod(i,subcol)==0)*subcol) = s;
        i=i+1;
        OutFitCellInner = {fit_cell{1,f}; RangedFitOut; FitTitles{f}};
        xlabel(AxesCell{3,2})
        ylabel(AxesCell{3,3})
        if ~isempty(p.Results.title)
            title([p.Results.title ' ' cd{3} ' - ' FitTitles{f}])
        else
            title([cd{3} ' - ' FitTitles{f}])
        end
        if ~isempty(p.Results.plt_cmds)
            for c = p.Results.plt_cmds(2,contains(p.Results.plt_cmds(1,:),cd{1},'IgnoreCase',true) | cellfun(@(c) isempty(c),p.Results.plt_cmds(1,:)))
                eval(c{:});
            end
        end
    end
    OutStruct.data = OutFitCellInner;
    OutStruct.axes = AxesCell(2,2:3);
    OutStruct.axesvec = axesvec;
    OutStruct.order = AxesCell(1,2:3);
    OutStruct.metadata.span = sprintf('%0.3g ',span);
    OutStruct.metadata.fitParams = StrFitUscore;
    OutStruct.metadata.fitNames = FitTitles;
    OutStructArray = [OutStructArray OutStruct];
end
OutAppendPath = [cd{4} '-' fit_axis '\' StrFitUscore ' ' p.Results.str StartTime];
set(fig,'Position',get(0,'Screensize'));
end

