function [sbp] = fit_cell_3D(data_cell, axes_cell, fit_func, fit_axis, varargin)
CharOrString = @(s) ischar(s) || isstring(s);
p = inputParser;
p.KeepUnmatched=true;
% p.addParameter('range', {}, @iscell);
p.addParameter('plt_select_data', {}, @iscell);
p.addParameter('plt_cmds', {}, @iscell);
p.addParameter('title', '', CharOrString);
p.addParameter('plot_fit', []);
p.parse(varargin{:});

d123 = 1:3;
ftAxPos = strcmpi(axes_cell(1,:),fit_axis);
perm_vec = [d123(ftAxPos) d123(~ftAxPos)];
axes_cell(2,:) = cellfun(@(c) permute(c,perm_vec), axes_cell(2,:) , 'UniformOutput', false);
axes_cell(:,:) = axes_cell(:,perm_vec);
data_cell(2,:) = cellfun(@(c) permute(c,perm_vec), data_cell(2,:) , 'UniformOutput', false);

FitAxMat = axes_cell{2,1};
AxJ = shiftdim(mean(axes_cell{2,2}));
AxK = shiftdim(mean(axes_cell{2,3}));
FitOut = zeros(size(FitAxMat, [2 3]));

if ~isempty(p.Results.plt_select_data)
    data_cell = data_cell(:,cellfun(@(c) contains(c,p.Results.plt_select_data(1,:)),data_cell(1,:)));
end

cd_mat = data_cell{2,1};
[~, ~, ~, fit_cell] = fit_func(FitAxMat(:,round(size(FitAxMat,2)/2),round(size(FitAxMat,3)/2)), cd_mat(:,round(size(FitAxMat,2)/2),round(size(FitAxMat,3)/2)));
fmax = size(fit_cell,2);
StrFitUscore = [];
for c = fit_cell(1,:)
    StrFitUscore = [StrFitUscore '_' c{:}];
end
StrFitUscore(1) = '';
FitTitles = fit_cell(3,:);
[subrow, subcol] = subplot_min_rectangle(size(data_cell,2)*fmax);
fig = figure;
i=1;

for cd = data_cell
    cd_mat = cd{2};
    for j=1:size(cd_mat, 2)
        for k=1:size(cd_mat, 3)
            [fit_str, OutFitFunc, span, fit_cell] = fit_func(FitAxMat(:,j,k), cd_mat(:,j,k));
            fit_progress = ['j-k = ' num2str(j) '-' num2str(k) ' ' axes_cell{4,2} ' = ' num2str(AxJ(j,1),3) ' ' axes_cell{4,3} ' = ' num2str(AxK(j,k),3)];
            if contains(fit_str, 'complex','IgnoreCase',true)
                disp(fit_progress)
            end
            if ~isempty(p.Results.plot_fit)
                if isfield(p.Results.plot_fit, 'func')
                    plot_func = p.Results.plot_fit.func;
                else
                    plot_func = @plot;
                end
                ax = s_plot(FitAxMat(:,j,k), cd_mat(:,j,k), '', [p.Results.title ' ' cd{4} ' ' strrep(fit_progress,'_','\_')], fit_str, '', axes_cell{3,1}, cd{3}, 10, plot_func);
                if ~isempty(OutFitFunc)
                    hold(ax(1),'on');
                    fplot(OutFitFunc, span)
                    hold(ax(1),'off');
                end
                if isfield(p.Results.plot_fit, 'savepath')
                    if ~isempty(OutFitFunc)
                        RealSavePath = [p.Results.plot_fit.savepath '\' cd{4} '\' StrFitUscore '\success' ];
                    else
                        RealSavePath = [p.Results.plot_fit.savepath '\' cd{4} '\' StrFitUscore '\fail' ];
                    end
                    if ~exist(RealSavePath, 'dir')
                        mkdir(RealSavePath)
                    end
                    saveas(ax(1), [RealSavePath '\' [p.Results.title ' ' cd{4} ' ' fit_progress] '.png']);
                end
            end
            for f=1:fmax
                FitOut(j,k,f) = fit_cell{2,f};
            end            
        end
    end
    for f=1:fmax
        figure(fig);
        s = subplot(subrow, subcol, i);
        sbp(subrow, subcol) = s;
        i=i+1;
        surf(AxJ, AxK, FitOut(:,:,f));
        xlabel(axes_cell{3,2})
        ylabel(axes_cell{3,3})
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
end          
end

