%% Load data file (struct 3D format) and get axes order
clear
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
DataStruct0 = LoadedFile.(fNames{contains(fNames, SetName)});
Order = DataStruct0.order;
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
[PlotFig, SubPlots, DataStruct] = process_plot_struct_data3D(DataStruct0, Order, slice_planes, 'title', title, 'ax_range', ax_range, 'val_range', val_range, 'plt_select_data', plt_select_data, 'plt_cmds', plt_cmds);

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
clear 'DataFigureHandles'
k = 0;
for j = 1:size(SubPlots,2)
    for i = 1:size(SubPlots,1)
        k = k+1;
        DataFigureHandles(k).fig = figure('Position',get(0,'Screensize'));
        DataFigureHandles(k).ax3D = subplot(1,3,1);
        DataFigureHandles(k).ax2D = subplot(1,3,2);
        DataFigureHandles(k).ax1D = subplot(1,3,3);
        
        copyobj(get(SubPlots(i,j), 'Children'),DataFigureHandles(k).ax3D)
        DataFigureHandles(k).ax3D = update_structure(DataFigureHandles(k).ax3D, SubPlots(i,j), 'ignore',{'Parent'}, 'new', true);
        DataFigureHandles(k).ax3D.Position = [0.05 DataFigureHandles(k).ax3D.Position(2:end)];
        colorbar(DataFigureHandles(k).ax3D, 'north')
        
        P3D = DataFigureHandles(k).ax3D.Position;
        Temp = DataFigureHandles(k).ax3D.UserData.axesvec;
        Tempf = fieldnames(Temp)';
        ax_lbls = {Temp.(Tempf{1}).label Temp.(Tempf{2}).label Temp.(Tempf{3}).label};
        ax_ranges = {[min(Temp.(Tempf{1}).vec) max(Temp.(Tempf{1}).vec)] [min(Temp.(Tempf{2}).vec) max(Temp.(Tempf{2}).vec)] [min(Temp.(Tempf{3}).vec) max(Temp.(Tempf{3}).vec)]};
        PopStrVal = [ax_lbls ; ax_ranges ; Tempf];
        PopOne2D = PopStrVal{2,2};
        PopOne3D = PopStrVal{2,1};
        
        P2D = DataFigureHandles(k).ax2D.Position;
        SliderCallFromPop2D = @(source,callbackdata,AxisString) slider_callback2D(source,callbackdata,DataFigureHandles(k).ax1D,AxisString,DataFigureHandles(k).ax2D);
        DataFigureHandles(k).s2D = uicontrol('Parent',DataFigureHandles(k).fig,'Style','slider', 'Units','normalized', 'Callback', @(src,cbdta) SliderCallFromPop2D(src,cbdta,PopStrVal{3,2}),...
            'value',mean(PopOne2D), 'min',PopOne2D(1), 'max',PopOne2D(2));
        DataFigureHandles(k).s2D.Position = [P2D(1) P2D(2)-0.08 P2D(3) 0.02];
        LeftTxt2D = uicontrol('Style','text', 'Units','normalized', 'Position',[P2D(1)-0.02 P2D(2)-0.08 0.02 0.02],'String',num2str(PopOne2D(1)) ,'FontSize',15);
        RightTxt2D = uicontrol('Style','text', 'Units','normalized', 'Position',[P2D(1)+P2D(3) P2D(2)-0.08 0.02 0.02],'String',num2str(PopOne2D(2)) ,'FontSize',15);
        
        PopCall2D = @(source,callbackdata,PopStrVal) popmenu_callback2D(source,callbackdata,DataFigureHandles(k).s2D,SliderCallFromPop2D,PopStrVal,LeftTxt2D,RightTxt2D);
        DataFigureHandles(k).p2D = uicontrol('Parent',DataFigureHandles(k).fig, 'Style','popupmenu', 'string',PopStrVal(1,2:3), 'Units','normalized', 'Position', [P2D(1)+0.12 P2D(2)-0.09 0.10 0.05], 'Callback',@(src,cbdta)PopCall2D(src,cbdta,PopStrVal(:,2:3)));
        
        SliderCallFromPop3D = @(source,callbackdata,AxisString) slider_callback3D(source,callbackdata,DataFigureHandles(k).ax2D,AxisString,DataFigureHandles(k).ax3D.UserData);
        DataFigureHandles(k).s3D = uicontrol('Parent',DataFigureHandles(k).fig,'Style','slider', 'Units','normalized', 'Callback', @(src,cbdta) SliderCallFromPop3D(src,cbdta,PopStrVal{3,1}),...
            'value',mean(PopOne3D), 'min',PopOne3D(1), 'max',PopOne3D(2));
        DataFigureHandles(k).s3D.Position = [P3D(1) P3D(2)-0.08 P3D(3) 0.02];
        LeftTxt3D = uicontrol('Style','text', 'Units','normalized', 'Position',[P3D(1)-0.02 P3D(2)-0.08 0.02 0.02],'String',num2str(PopOne3D(1)) ,'FontSize',15);
        RightTxt3D = uicontrol('Style','text', 'Units','normalized', 'Position',[P3D(1)+P3D(3) P3D(2)-0.08 0.02 0.02],'String',num2str(PopOne3D(2)) ,'FontSize',15);
        
        PopCall3D = @(source,callbackdata) popmenu_callback3D(source,callbackdata,DataFigureHandles(k).s3D,SliderCallFromPop2D,DataFigureHandles(k).p2D,PopStrVal,LeftTxt3D,RightTxt3D);
        DataFigureHandles(k).p3D = uicontrol('Parent',DataFigureHandles(k).fig, 'Style','popupmenu', 'string',PopStrVal(1,:), 'Units','normalized', 'Position', [P3D(1)+0.12 P3D(2)-0.09 0.10 0.05], 'Callback',PopCall3D);
        
        
        
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
CV_plot_fit.properties = {'Marker','o', 'LineStyle','none'};
% CV_plot_fit = [];

CV_fit_select = {'capacitance'};
CV_fit_axis = 'offset';

CV_val_range = {};
CV_val_range(:,end+1) = {''; [0 inf]};
CV_val_range(:,end+1) = {''; [0 inf]};
%Fit function and options
CV_offset_range = [-inf inf];
CV_str = ['fit rng=' num2str(CV_offset_range) ' '];
syms V N Vb n
limits_A = {'Vb',[0 0.5 1],'N',[1e16 2e18 5e18],'n',[1 3 10]};
IdealityFactor = 1;
es = 11.68;
A = (150^2 * pi)*1e-8;
CV_fit_func = @(x,y) C_schot_fit_A(x,y,CV_offset_range,A,es,N,Vb,IdealityFactor);

[CV_fig, CV_sbp, CV_struct_array, CV_FitAppendPath, CV_FitName] = fit_cell_3D(DataStruct.data, DataStruct.axes, CV_fit_func, CV_fit_axis, 'plt_select_data',CV_fit_select, 'title',title, 'plot_fit',CV_plot_fit, 'val_range',CV_val_range, 'str',CV_str);

if Save_CV_plot
    CV_save_path = [SavePlotPath '\' CV_FitAppendPath];
    mkdir(CV_save_path);
    save([CV_save_path '\CV fit out ' CV_FitName], 'CV_struct_array')
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
DLCP_plot_fit.properties = {'Marker','o', 'LineStyle','none'};
DLCP_plot_fit = [];

DLCP_fit_select = {'capacitance'};  
DLCP_fit_axis = 'amplitude'; 

DLCP_val_range = {};
DLCP_val_range(:,end+1) = {'N'; [0 1e18]};
%Fit function and options
DLCP_amplitude_Range = [0.05 inf];
DLCPstr = ['fit rng=' num2str(DLCP_offset_range) ' '];
A = (150^2 * pi)*1e-8;
DLCP_fit_func = @(x,y) DLCP_fit(x,y, A, 'range',DLCP_amplitude_Range);

[DLCP_fig, DLCP_sbp, DLCP_struct_array, DLCP_FitAppendPath, DLCP_FitName] = fit_cell_3D(DataStruct.data, DataStruct.axes, DLCP_fit_func, DLCP_fit_axis, 'plt_select_data',DLCP_fit_select, 'title',title, 'plot_fit',DLCP_plot_fit, 'val_range', DLCP_val_range, 'str',DLCP_str);
if Save_DLCP_plot
    DLCP_save_path = [SavePlotPath '\' DLCP_FitAppendPath];
    mkdir(DLCP_save_path);
    save([DLCP_save_path '\DLCP fit_ ut ' DLCP_FitName], 'DLCP_struct_array')
    f=figure('Position',get(0,'Screensize'));
    copyobj(get(DLCP_fig, 'Children'),f)
    f.Visible = 'on';
    saveas(f, [DLCP_save_path '\DLCP fit ' DLCP_FitName '.fig']);
    f.Visible = 'off';
    saveas(f, [DLCP_save_path '\DLCP fit ' DLCP_FitName '.bmp']);
    close(f)
end

%% Functions

function order = get_order_from_path(path, varargin)
if isempty(varargin)
    search_order = 'order_'; % input the string just before the first parameter
else
    search_order = varargin{:};
end
order_string = path((strfind(path,search_order)+length(search_order)):end);
order_ind = strfind(order_string, '_');
order = {order_string(1:order_ind(1)-1), order_string(order_ind(1)+1:order_ind(2)-1), order_string(order_ind(2)+1:end)};
end

function slider_callback3D(src,cbdata,Plot2D,SelectAxString,DataStruct)
    Axes = DataStruct.axesvec;
    AxNames = fieldnames(Axes);
    SelectAxN = find(strcmpi(AxNames,SelectAxString));
    AxNames = AxNames(~strcmpi(AxNames,SelectAxString));
    [~,SelectMinInd] = min(abs(Axes.(SelectAxString).vec - src.Value));
    switch SelectAxN
        case 1
            C = DataStruct.data(SelectMinInd,:,:);
        case 2
            C = DataStruct.data(:,SelectMinInd,:);
        case 3
            C = DataStruct.data(:,:,SelectMinInd);
    end
    if ndims(C)==3 
        if SelectAxN==2
            C = shiftdim(C,2)';
        else
            C = shiftdim(C);
        end
    end
        
    pcolor(Axes.(AxNames{2}).vec, Axes.(AxNames{1}).vec, C, 'Parent',Plot2D);
    colormap(Plot2D, 'turbo')
    colorbar(Plot2D, 'northoutside')
    Plot2D.XLabel.String = Axes.(AxNames{2}).label;
    Plot2D.YLabel.String = Axes.(AxNames{1}).label;
    Plot2D.Title.String =[Axes.(SelectAxString).label ' = ' num2str(Axes.(SelectAxString).vec(SelectMinInd))];
    
end

function slider_callback2D(src,cbdata,Plot1D,AxisString,Plot2D)
    P2DC = Plot2D.Children;
    P2DX = Plot2D.XLabel.String;
    P2DY = Plot2D.YLabel.String;
    if strcmpi(P2DX,AxisString)
        [~,SelectMinInd] = min(abs(P2DC.XData - src.Value));
        plot(P2DC.YData, shiftdim(P2DC.CData(:,SelectMinInd)))
        Plot1D.Title.String =[P2DX ' = ' num2str(P2DC.XData(SelectMinInd))];
        Plot1D.XLabel.String = P2DY;
    elseif strcmpi(P2DY,AxisString)
        [~,SelectMinInd] = min(abs(P2DC.YData - src.Value));
        plot(P2DC.XData, shiftdim(P2DC.CData(SelectMinInd,:)))
        Plot1D.Title.String =[P2DY ' = ' num2str(P2DC.YData(SelectMinInd))];
        Plot1D.XLabel.String = P2DX;
    end 
    
end

function popmenu_callback3D(src,cbdata,slider3D,SliderCall,Pop2D,Pop2D_call,value_cell,LeftTxt,RightTxt)
    SV = src.Value;
    AxBounds = value_cell{2,SV};
    set(slider3D, 'value',mean(AxBounds), 'min',AxBounds(1), 'max',AxBounds(2), 'Callback', @(src,cbdata) SliderCall(src,cbdata,value_cell{3,SV}))
    set(LeftTxt, 'String', num2str(AxBounds(1)))
    set(RightTxt, 'String', num2str(AxBounds(2)))
    value_cell(:,SV) = [];
    set(Pop2D, 'string',value_cell(1,:), 'Callback',@(src,cbdata)Pop2D_call(src,cbdata,value_cell));
end

function popmenu_callback2D(src,cbdata,slider2D,SliderCall,value_cell,LeftTxt,RightTxt)
    SV = src.Value;
    AxBounds = value_cell{2,SV};
    set(slider2D, 'value',mean(AxBounds), 'min',AxBounds(1), 'max',AxBounds(2), 'Callback', @(src,cbdata) SliderCall(src,cbdata,value_cell{1,SV}))
    set(LeftTxt, 'String', num2str(AxBounds(1)))
    set(RightTxt, 'String', num2str(AxBounds(2)))
end

