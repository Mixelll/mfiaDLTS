function DFH = MFIA_3D_PGUI(SubPlots)
clear 'DataFigureHandles'
k = 0;
Nsbp = numel(SubPlots);
for j = 1:size(SubPlots,2)
    for i = 1:size(SubPlots,1)
        k = k+1;
        DFH(k).fig = figure('Position',get(0,'Screensize'));
        DFH(k).ax3D = subplot(1,3,1);
        DFH(k).ax2D = subplot(1,3,2); 
        DFH(k).ax1D = subplot(1,3,3);
        
        copyobj(get(SubPlots(i,j), 'Children'),DFH(k).ax3D)
        DFH(k).ax3D = update_structure(DFH(k).ax3D, SubPlots(i,j), 'ignore',{'Parent'}, 'new', true);
        DFH(k).ax2D.UserData.type = DFH(k).ax3D.UserData.type;
        DFH(k).ax1D.UserData.type = DFH(k).ax3D.UserData.type;
        if Nsbp>2 || Nsbp==1
            DFH(k).ax3D.Position = [0.05    0.11    0.304    0.815];
        else
            DFH(k).ax3D.Position = [0.05 DFH(k).ax3D.Position(2:end)];
        end
        colorbar(DFH(k).ax3D, 'north')
        
        P3D = DFH(k).ax3D.Position;
        Temp = DFH(k).ax3D.UserData.axesvec;
        Tempf = fieldnames(Temp)';
        ax_lbls = {Temp.(Tempf{1}).label Temp.(Tempf{2}).label Temp.(Tempf{3}).label};
        ax_ranges = {[min(Temp.(Tempf{1}).vec) max(Temp.(Tempf{1}).vec)] [min(Temp.(Tempf{2}).vec) max(Temp.(Tempf{2}).vec)] [min(Temp.(Tempf{3}).vec) max(Temp.(Tempf{3}).vec)]};
        ax_lengths = {length(Temp.(Tempf{1}).vec) length(Temp.(Tempf{2}).vec) length(Temp.(Tempf{3}).vec)};
        PopStrVal = [ax_lbls ; ax_ranges ; Tempf; ax_lengths];
        PopOne2D = PopStrVal{2,2};
        PopOne3D = PopStrVal{2,1};
        SliderStep2 = 0.1;
        
        P2D = DFH(k).ax2D.Position;
        SliderCallFromPop2D = @(source,callbackdata,AxisString) slider_callback2D(source,callbackdata,DFH(k).ax1D,AxisString,DFH(k).ax2D);
        DFH(k).s2D = uicontrol('Parent',DFH(k).fig,'Style','slider', 'Units','normalized', 'Callback', @(src,cbdta) SliderCallFromPop2D(src,cbdta,PopStrVal{3,2}),...
            'value',mean(PopOne2D), 'min',PopOne2D(1), 'max',PopOne2D(2), 'SliderStep',[1/PopStrVal{4,2} SliderStep2]);
        DFH(k).s2D.Position = [P2D(1) P2D(2)-0.08 P2D(3) 0.02];
        LeftTxt2D = uicontrol('Style','text', 'Units','normalized', 'Position',[P2D(1)-0.02 P2D(2)-0.08 0.02 0.02],'String',num2str(PopOne2D(1)) ,'FontSize',15);
        RightTxt2D = uicontrol('Style','text', 'Units','normalized', 'Position',[P2D(1)+P2D(3) P2D(2)-0.08 0.02 0.02],'String',num2str(PopOne2D(2)) ,'FontSize',15);
        
        PopCall2D = @(source,callbackdata,PopStrVal) popmenu_callback2D(source,callbackdata,DFH(k).s2D,SliderCallFromPop2D,PopStrVal,LeftTxt2D,RightTxt2D);
        DFH(k).p2D = uicontrol('Parent',DFH(k).fig, 'Style','popupmenu', 'string',PopStrVal(1,2:3), 'Units','normalized', 'Position', [P2D(1)+0.08 P2D(2)-0.09 0.08 0.05], 'Callback',@(src,cbdta)PopCall2D(src,cbdta,PopStrVal(:,2:3)));
        
        SliderCallFromPop3D = @(source,callbackdata,AxisString) slider_callback3D(source,callbackdata,DFH(k).ax2D,AxisString,DFH(k).ax3D.UserData);
        DFH(k).s3D = uicontrol('Parent',DFH(k).fig,'Style','slider', 'Units','normalized', 'Callback', @(src,cbdta) SliderCallFromPop3D(src,cbdta,PopStrVal{3,1}),...
            'value',mean(PopOne3D), 'min',PopOne3D(1), 'max',PopOne3D(2), 'SliderStep',[1/PopStrVal{4,1} SliderStep2]);
        DFH(k).s3D.Position = [P3D(1) P3D(2)-0.08 P3D(3) 0.02];
        LeftTxt3D = uicontrol('Style','text', 'Units','normalized', 'Position',[P3D(1)-0.02 P3D(2)-0.08 0.02 0.02],'String',num2str(PopOne3D(1)) ,'FontSize',15);
        RightTxt3D = uicontrol('Style','text', 'Units','normalized', 'Position',[P3D(1)+P3D(3) P3D(2)-0.08 0.02 0.02],'String',num2str(PopOne3D(2)) ,'FontSize',15);
        
        PopCall3D = @(source,callbackdata) popmenu_callback3D(source,callbackdata,DFH(k).s3D,SliderCallFromPop3D,DFH(k).p2D,PopCall2D,PopStrVal,LeftTxt3D,RightTxt3D);
        DFH(k).p3D = uicontrol('Parent',DFH(k).fig, 'Style','popupmenu', 'string',PopStrVal(1,:), 'Units','normalized', 'Position', [P3D(1)+0.11 P3D(2)-0.09 0.08 0.05], 'Callback',PopCall3D);
         
    end
end

%% Functions
function slider_callback3D(src,cbdata,Plot2D,SelectAxString,DataStruct)
    UserData = Plot2D.UserData;
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
    Plot2D.UserData = UserData;
end

function slider_callback2D(src,cbdata,Plot1D,AxisString,Plot2D)
    UserData = Plot1D.UserData;
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
    Plot1D.UserData = UserData;
    Plot1D.YLabel.String = Plot1D.UserData.type;
end

function popmenu_callback3D(src,cbdata,slider3D,SliderCall,Pop2D,Pop2D_call,value_cell,LeftTxt,RightTxt)
    SV = src.Value;
    AxBounds = value_cell{2,SV};
    set(slider3D, 'value',mean(AxBounds), 'min',AxBounds(1), 'max',AxBounds(2), 'SliderStep',[1/value_cell{4,SV} slider3D.SliderStep(2)],...
        'Callback', @(src,cbdata) SliderCall(src,cbdata,value_cell{3,SV}))
    set(LeftTxt, 'String', num2str(AxBounds(1)))
    set(RightTxt, 'String', num2str(AxBounds(2)))
    value_cell(:,SV) = [];
    set(Pop2D, 'string',value_cell(1,:), 'Callback',@(src,cbdata)Pop2D_call(src,cbdata,value_cell));
end

function popmenu_callback2D(src,cbdata,slider2D,SliderCall,value_cell,LeftTxt,RightTxt)
    SV = src.Value;
    AxBounds = value_cell{2,SV};
    set(slider2D, 'value',mean(AxBounds), 'min',AxBounds(1), 'max',AxBounds(2),  'SliderStep',[1/value_cell{4,SV} slider2D.SliderStep(2)],...
        'Callback', @(src,cbdata) SliderCall(src,cbdata,value_cell{1,SV}))
    set(LeftTxt, 'String', num2str(AxBounds(1)))
    set(RightTxt, 'String', num2str(AxBounds(2)))
end
end

