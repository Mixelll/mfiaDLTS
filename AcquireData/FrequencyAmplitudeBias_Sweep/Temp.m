%% Load data file (struct 3D format) and get axes order
FolderPath = ['C:\Users\' getenv('USERNAME') '\Google Drive\MATLAB Drive\vars'];
% FolderPath = 'D:\Google Drive\MATLAB Drive\vars';
SampleName = 'B5 b5 150um 9';
SweepID = 'A';
FileName = [SampleName ' ' SweepID];
SetName = '3D_sweep_order';
RealPath = [FolderPath '\' SampleName '\' FileName];
%RealPath = '' % Or put the path yourself

LoadedFile = load(RealPath);
fNames = fieldnames(LoadedFile);
DataStruct = LoadedFile.(fNames{contains(fNames, SetName)});
Order = DataStruct.order;
% StrBeforeOrder = 'order_'; % input the string just before the first parameter
% Order = MFIA_get_order_from_path(fNames{contains(fNames,SetName)}, StrBeforeOrder)
title = RealPath(sum(find(RealPath=='\', 1, 'last'))+1:end);
SavePlotPath = [FolderPath '\' SampleName '\' title];


%% Plot and Save data
%Save data slice plots
SaveDataPlots = false;
%Select data to plot (leave empty to plot all data in 3D struct) and define plot aliases (from IA parametric model)
plt_select_data = {};
plt_select_data(:,end+1) = {'param0'; 'Resistance [ohm]'};
plt_select_data(:,end+1) = {'param1'; 'Capacitance [F]'};
%Select axis range - data outside range is discarded from here onward (applies to all subsequent fits)
ax_range = {};
ax_range(:,end+1) = {'frequency'; [0 inf]};
ax_range(:,end+1) = {'amplitude'; [0 inf]}; 
ax_range(:,end+1) = {'offset'; [-inf inf]};
%Select value range - data outside range is replaced with NaN from here onward (applies to all subsequent fits)
val_range = {};
val_range(:,end+1) = {'param0'; [0 inf]};
val_range(:,end+1) = {'param1'; [0 inf]};
%Define slice planes
slice_planes = {};
slice_planes(:,end+1) = {'frequency'; [1e2 1e3 1e4 1e5 5e5]};
slice_planes(:,end+1) = {'amplitude'; []}; 
slice_planes(:,end+1) = {'offset'; [0 -0.1 -0.2 -0.5 -1 -1.5]};
%Add plot formatting commands: All: plt_cmds(2,end+1) =  {'comand'}; '' for ' inside string. Target: plt_cmds(:,end+1) =  {'target' ; 'comand'};
plt_cmds = {};
% plt_cmds(2,end+1) = {'grid on'};
plt_cmds(2,end+1) = {'colorbar(''eastoutside'')'};
plt_cmds(2,end+1) = {'colormap(s,''turbo'')'};
plt_cmds(2,end+1) = {'colormap(s,interp1(colormap(s), 1:(length(colormap(s))/998):length(colormap(s))))'};
plt_cmds(:,end+1) = {'param0'; 'colormap(s,[compress_array_exp(colormap(s),10,10, ''Reverse'', true, ''Interp'', true)])'};
% plt_cmds(2,end+1) = {'colormap(s,[0 0 0; colormap(s)])'};
% plt_cmds(2,end+1) = {'colormap([colormap; 1 1 1])'};
plt_cmds(:,end+1) = {'param0'; 'caxis([0 1e7])'};
plt_cmds(:,end+1) = {'param1'; 'caxis([0 inf])'};

plt_log_freq = true; % true for log plot frequency
if plt_log_freq
    Order{contains(Order, 'frequency')} = 'log_frequency';
    ax_range(:,contains(ax_range(1,:), 'frequency')) = {'log_frequency'; log10(ax_range{2,contains(ax_range(1,:), 'frequency')})};
    slice_planes(:,contains(slice_planes(1,:), 'frequency')) = {'log_frequency'; log10(slice_planes{2,contains(slice_planes(1,:), 'frequency')})};
end
%Plot slice planes of 3D data
[PlotFig, SubPlots, DataCell, AxesCell] = process_plot_struct_data3D(DataStruct, Order, slice_planes, 'title', title, 'ax_range', ax_range, 'val_range', val_range, 'plt_select_data', plt_select_data, 'plt_cmds', plt_cmds);

if SaveDataPlots
    StartTime = datestr(now, 'yyyy-mm-dd HH-MM-SS');
    mkdir(SavePlotPath);
    f=figure('Position',get(0,'Screensize'));
    copyobj(get(PlotFig, 'Children'),f)
    f.Visible = 'on';
    saveas(f, [SavePlotPath '\' FileName '_slice_plot ' StartTime '.fig']);
    f.Visible = 'off';
    saveas(f, [SavePlotPath '\' FileName '_slice_plot ' StartTime '.bmp']);
    close(f)
end
%         DataSubPlotHandles=figure('Position',get(0,'Screensize'));
clear 'DataFigureHandles(k)'
k = 0;
for j = 1:size(SubPlots,2)
    for i = 1:size(SubPlots,1)
        k = k+1;
        DataFigureHandles(k).fig = figure;
        DataFigureHandles(k).ax3D = subplot(1,3,1);
        copyobj(get(SubPlots(i,j), 'Children'),DataFigureHandles(k).ax3D)
        DataFigureHandles(k).ax3D = update_structure(DataFigureHandles(k).ax3D, SubPlots(i,j), 'ignore',{'Parent'}, 'new', false);
        DataFigureHandles(k).ax3D.Position = [0.05 DataFigureHandles(k).ax3D.Position(2:end)];
        P3D = DataFigureHandles(k).ax3D.Position;
        PopStrVal = {DataFigureHandles(k).ax3D.XLabel.String, DataFigureHandles(k).ax3D.YLabel.String, DataFigureHandles(k).ax3D.ZLabel.String; DataFigureHandles(k).ax3D.XLim, DataFigureHandles(k).ax3D.YLim, DataFigureHandles(k).ax3D.ZLim};
        
        DataFigureHandles(k).p3D = uicontrol('Parent',DataFigureHandles(k).fig, 'Style','popupmenu', 'string',PopStrVal(1,:), 'Units','normalized', 'Position', [P3D(1) P3D(2)-0.1 P3D(3) 0.05]);
        
        DataFigureHandles(k).s3D = uicontrol('Parent',DataFigureHandles(k).fig,'Style','slider','Position',[P3D(1) P3D(2)-0.1 P3D(3) 0.05],...
               'Units','normalized', 'value',mean(PopStrVal{2,DataFigureHandles(k).p3D.Value}), 'min',PopStrVal{2,DataFigureHandles(k).p3D.Value}(1), 'max',PopStrVal{2,DataFigureHandles(k).p3D.Value}(2));
        DataFigureHandles(k).s3D.callback = @(es,ed) updateSystem(h,DataFigureHandles(k).s3D.Value); 
        DataFigureHandles(k).ax2D = subplot(1,3,2);
        DataFigureHandles(k).ax1D = subplot(1,3,3);
    end
end

%% Fit data along a certain dimension
%C-V fit
%Plot and save options
Save_CV_plot = true;

CV_plot_fit = [];
CV_plot_fit.func = @plot;
CV_plot_fit.visible = 0;
CV_plot_fit.savepath = SavePlotPath;
% CV_plot_fit = [];

CV_fit_select = {'capacitance'};
CV_fit_axis = 'offset';

CV_val_range = {};
CV_val_range(:,end+1) = {''; [0 inf]};
CV_val_range(:,end+1) = {''; [0 inf]};
%Fit function and options
CV_offset_range = [-inf inf];
syms V N Vb n
limits_A = {'Vb',[0 0.5 1],'N',[1e16 2e18 5e18],'n',[1 3 10]};
IdealityFactor = 1;
es = 11.68;
A = (150^2 * pi)*1e-8;
CV_fit_func = @(x,y) C_schot_fit_A(x,y,CV_offset_range,A,es,N,Vb,IdealityFactor);

[CV_fig, CV_sbp, CV_struct_array, CV_data_cell, CV_axes_cell, CV_FitAppendPath, CV_FitName] = fit_cell_3D(DataCell, AxesCell, CV_fit_func, CV_fit_axis, 'plt_select_data',CV_fit_select, 'title',title, 'plot_fit',CV_plot_fit, 'val_range', CV_val_range);

if Save_CV_plot
    CV_save_path = [SavePlotPath '\' CV_FitAppendPath];
    mkdir(CV_save_path);
    save([CV_save_path '\CV fit out ' CV_FitName], 'CV_struct_array','CV_data_cell','CV_axes_cell')
    f=figure('Position',get(0,'Screensize'));
    copyobj(get(CV_fig, 'Children'),f)
    f.Visible = 'on';
    saveas(f, [CV_save_path '\CV fit ' CV_FitName '.fig']);
    f.Visible = 'off';
    saveas(f, [CV_save_path '\CV fit ' CV_FitName '.bmp']);
    close(f)
end

%% DLCP fit
Save_DLCP_plot = true;

DLCP_plot_fit = [];
DLCP_plot_fit.func = @plot;
DLCP_plot_fit.visible = 0;
DLCP_plot_fit.savepath = SavePlotPath;
DLCP_plot_fit = [];

DLCP_fit_select = {'capacitance'};  
DLCP_fit_axis = 'amplitude'; 

DLCP_val_range = {};
DLCP_val_range(:,end+1) = {'N'; [0 1e18]};
%Fit function and options
DLCP_amplitude_Range = [0.05 inf];
A = (150^2 * pi)*1e-8;
DLCP_fit_func = @(x,y) DLCP_fit(x,y, A, 'range',DLCP_amplitude_Range);

[DLCP_fig, DLCP_sbp, DLCP_struct_array, DLCP_data_cell, DLCP_axes_cell, DLCP_FitAppendPath, DLCP_FitName] = fit_cell_3D(DataCell, AxesCell, DLCP_fit_func, DLCP_fit_axis, 'plt_select_data',DLCP_fit_select, 'title',title, 'plot_fit',DLCP_plot_fit, 'val_range', DLCP_val_range);
if Save_DLCP_plot
    DLCP_save_path = [SavePlotPath '\' DLCP_FitAppendPath];
    mkdir(DLCP_save_path);
    save([DLCP_save_path '\DLCP fit_ ut ' DLCP_FitName], 'DLCP_struct_array','DLCP_data_cell','DLCP_axes_cell')
    f=figure('Position',get(0,'Screensize'));
    copyobj(get(DLCP_fig, 'Children'),f)
    f.Visible = 'on';
    saveas(f, [DLCP_save_path '\DLCP fit ' DLCP_FitName '.fig']);
    f.Visible = 'off';
    saveas(f, [DLCP_save_path '\DLCP fit ' DLCP_FitName '.bmp']);
    close(f)
end
