function [figs] = MFIA_PGUI(SubPlots,D)
clear 'DataFigureHandles'
k = 0;
Nsbp = numel(SubPlots);
if D ==3
    for j = 1:size(SubPlots,2)
        for i = 1:size(SubPlots,1)
            k = k+1;
            fig = figure('Position',get(0,'Screensize'), 'UserData',SubPlots(i,j).Title.UserData);
            ax3D = subplot(1,3,1);
            ax2D = subplot(1,3,2); 
            ax1D = subplot(1,3,3);

            copyobj(get(SubPlots(i,j), 'Children'),ax3D)
            ax3D = update_structure(ax3D, SubPlots(i,j), 'ignore',{'Parent'}, 'new', true);
            ax2D.UserData.type = ax3D.UserData.type;
            ax1D.UserData.type = ax3D.UserData.type;
            if Nsbp>2 || Nsbp==1
                ax3D.Position = [0.05    0.11    0.304    0.815];
            else
                ax3D.Position = [0.05 ax3D.Position(2:end)];
            end
            colorbar(ax3D, 'north')
            ax3D.Units = 'pixels'; ax2D.Units = 'pixels'; ax1D.Units = 'pixels';
            P3D = ax3D.Position; P2D = ax2D.Position; P1D = ax1D.Position;
            ax3D.Units = 'normalized'; ax2D.Units = 'normalized'; ax1D.Units = 'normalized';
            Temp = ax3D.UserData.axesvec;
            Tempf = fieldnames(Temp)';
            ax_lbls = {Temp.(Tempf{1}).label Temp.(Tempf{2}).label Temp.(Tempf{3}).label};
            ax_ranges = {[min(Temp.(Tempf{1}).vec) max(Temp.(Tempf{1}).vec)] [min(Temp.(Tempf{2}).vec) max(Temp.(Tempf{2}).vec)] [min(Temp.(Tempf{3}).vec) max(Temp.(Tempf{3}).vec)]};
            ax_lengths = {length(Temp.(Tempf{1}).vec) length(Temp.(Tempf{2}).vec) length(Temp.(Tempf{3}).vec)};
            PopStrVal = [ax_lbls ; ax_ranges ; Tempf; ax_lengths];
            PopOne2D = PopStrVal{2,2};
            PopOne3D = PopStrVal{2,1};
            SliderStep2 = 0.1;

            
            
            c1Dpc = uicontrol('Parent',fig,'Style','checkbox', 'Units','pixels', 'value',0, 'max',1, 'min',0, 'String','Keep Axes Lim', 'FontSize',13, 'Position',[P1D(1) P1D(2)-85 150 20]);
            c2Dpc = uicontrol('Parent',fig,'Style','checkbox', 'Units','pixels', 'value',0, 'max',1, 'min',0, 'String','Keep Axes Lim', 'FontSize',13, 'Position',[P2D(1) P2D(2)-85 150 20]);
            
            SliderCallFromPop2D = @(source,callbackdata,AxisString) slider_callback2D(source,callbackdata,ax1D,AxisString,ax2D,c1Dpc);
            s2D = uicontrol('Parent',fig,'Style','slider', 'Units','pixels', 'Callback', @(src,cbdta) SliderCallFromPop2D(src,cbdta,PopStrVal{3,2}),...
                'value',mean(PopOne2D), 'min',PopOne2D(1), 'max',PopOne2D(2), 'SliderStep',[1/PopStrVal{4,2} SliderStep2]);
            s2D.Position = [P2D(1) P2D(2)-120 P2D(3) 30];
            LeftTxt2D = uicontrol('Style','text', 'Units','pixels', 'Position',[P2D(1)-20 P2D(2)-120 70 30],'String',num2str(PopOne2D(1)) ,'FontSize',15);
            RightTxt2D = uicontrol('Style','text', 'Units','pixels', 'Position',[P2D(1)+P2D(3) P2D(2)-120 70 30],'String',num2str(PopOne2D(2)) ,'FontSize',15);

            PopCall2D = @(source,callbackdata,PopStrVal) popmenu_callback2D(source,callbackdata,s2D,SliderCallFromPop2D,PopStrVal,LeftTxt2D,RightTxt2D);
            p2D = uicontrol('Parent',fig, 'Style','popupmenu', 'string',PopStrVal(1,2:3), 'FontSize',10, 'Units','pixels', 'Position', [P2D(1)+150 P2D(2)-110 250 50], 'Callback',@(src,cbdta)PopCall2D(src,cbdta,PopStrVal(:,2:3)));

            SliderCallFromPop3D = @(source,callbackdata,AxisString) slider_callback3D(source,callbackdata,ax2D,AxisString,ax3D.UserData,c2Dpc);
            s3D = uicontrol('Parent',fig,'Style','slider', 'Units','pixels', 'Callback', @(src,cbdta) SliderCallFromPop3D(src,cbdta,PopStrVal{3,1}),...
                'value',mean(PopOne3D), 'min',PopOne3D(1), 'max',PopOne3D(2), 'SliderStep',[1/PopStrVal{4,1} SliderStep2]);
            s3D.Position = [P3D(1) P3D(2)-120 P3D(3) 30];
            LeftTxt3D = uicontrol('Style','text', 'Units','pixels', 'Position',[P3D(1)-20 P3D(2)-120 70 30],'String',num2str(PopOne3D(1)) ,'FontSize',15);
            RightTxt3D = uicontrol('Style','text', 'Units','pixels', 'Position',[P3D(1)+P3D(3) P3D(2)-120 70 30],'String',num2str(PopOne3D(2)) ,'FontSize',15);

            PopCall3D = @(source,callbackdata) popmenu_callback3D(source,callbackdata,s3D,SliderCallFromPop3D,p2D,PopCall2D,PopStrVal,LeftTxt3D,RightTxt3D);
            p3D = uicontrol('Parent',fig, 'Style','popupmenu', 'string',PopStrVal(1,:), 'FontSize',10, 'Units','pixels', 'Position', [P3D(1)+300 P3D(2)-113 250 50], 'Callback',PopCall3D);
            figs(k) = fig;
        end
    end
elseif D==2
    for j = 1:size(SubPlots,2)
        for i = 1:size(SubPlots,1)
            k = k+1;
            fig = figure('Position',get(0,'Screensize'), 'UserData',SubPlots(i,j).Title.UserData);
            ax2D = subplot(1,3,1);
            ax2Dpc = subplot(1,3,2); 
            ax1D = subplot(1,3,3);

            copyobj(get(SubPlots(i,j), 'Children'),ax2D)
            ax2D = update_structure(ax2D, SubPlots(i,j), 'ignore',{'Parent'}, 'new', true);
            try 
                pcolor(mean(ax2D.Children.XData,2), mean(ax2D.Children.YData,1), ax2D.Children.CData.', 'Parent',ax2Dpc);
            catch
                pcolor(mean(ax2D.Children.XData,2), mean(ax2D.Children.YData,1), ax2D.Children.CData, 'Parent',ax2Dpc);
            end
            ax2Dpc.XLabel.String = ax2D.XLabel.String;
            ax2Dpc.YLabel.String = ax2D.YLabel.String;
            ax2Dpc.UserData.type = ax2D.UserData.type;
            ax1D.UserData.type = ax2D.UserData.type;
            
            if Nsbp>2 || Nsbp==1
                ax2D.Position = [0.05    0.11    0.304    0.815];
            else
                ax2D.Position = [0.05 ax2D.Position(2:end)];
            end
            colorbar(ax2Dpc, 'northoutside')
            ax2D.Units = 'pixels'; ax2Dpc.Units = 'pixels'; ax1D.Units = 'pixels';
            P2D = ax2D.Position; P2Dpc = ax2Dpc.Position; P1D = ax1D.Position;
            ax2D.Units = 'normalized'; ax2Dpc.Units = 'normalized'; ax1D.Units = 'normalized';
            Temp = ax2D.UserData.axesvec;
            Tempf = fieldnames(Temp)';
            ax_lbls = {Temp.(Tempf{1}).label Temp.(Tempf{2}).label};
            ax_ranges = {[min(Temp.(Tempf{1}).vec) max(Temp.(Tempf{1}).vec)] [min(Temp.(Tempf{2}).vec) max(Temp.(Tempf{2}).vec)]};
            ax_lengths = {length(Temp.(Tempf{1}).vec) length(Temp.(Tempf{2}).vec)};
            PopStrVal = [ax_lbls ; ax_ranges ; Tempf; ax_lengths];
            PopOne2D = PopStrVal{2,1};
            SliderStep2 = 0.1;
            
            c1Dpc = uicontrol('Parent',fig,'Style','checkbox', 'Units','pixels', 'value',0, 'max',1, 'min',0, 'String','Keep Axes Lim', 'FontSize',13, 'Position',[P1D(1) P1D(2)-85 150 20]);
            
            SliderCallFromPop2D = @(source,callbackdata,AxisString) slider_callback2D(source,callbackdata,ax1D,AxisString,ax2Dpc,c1Dpc);
            s2D = uicontrol('Parent',fig,'Style','slider', 'Units','pixels', 'Callback', @(src,cbdta) SliderCallFromPop2D(src,cbdta,PopStrVal{3,2}),...
                'value',mean(PopOne2D), 'min',PopOne2D(1), 'max',PopOne2D(2), 'SliderStep',[1/PopStrVal{4,2} SliderStep2]);
            s2D.Position = [P2Dpc(1) P2Dpc(2)-120 P2Dpc(3) 30];
            LeftTxt2D = uicontrol('Style','text', 'Units','pixels', 'Position',[P2Dpc(1)+20 P2Dpc(2)-120 70 30],'String',num2str(PopOne2D(1)) ,'FontSize',15);
            RightTxt2D = uicontrol('Style','text', 'Units','pixels', 'Position',[P2Dpc(1)+P2Dpc(3) P2Dpc(2)-120 70 30],'String',num2str(PopOne2D(2)) ,'FontSize',15);

            PopCall2D = @(source,callbackdata,PopStrVal) popmenu_callback2D(source,callbackdata,s2D,SliderCallFromPop2D,PopStrVal,LeftTxt2D,RightTxt2D);
            p2D = uicontrol('Parent',fig, 'Style','popupmenu', 'string',PopStrVal(1,:), 'FontSize',10, 'Units','pixels', 'Position', [P2Dpc(1)+150 P2Dpc(2)-110 250 50], 'Callback',@(src,cbdta)PopCall2D(src,cbdta,PopStrVal(:,:)));

            figs(k) = fig;
        end
    end
end
%% Functions
function slider_callback3D(src,cbdata,Plot2D,SelectAxString,DataStruct,CheckBox)
    XlblBefore = Plot2D.XLabel.String;
    YlblBefore = Plot2D.YLabel.String;
    XLim = Plot2D.XLim;
    YLim = Plot2D.YLim;
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
    if strcmpi(Plot2D.XLabel.String,XlblBefore) && strcmpi(Plot2D.YLabel.String,YlblBefore) && CheckBox.Value
        Plot2D.XLim = XLim;
        Plot2D.YLim = YLim;
    end  
end

function slider_callback2D(src,cbdata,Plot1D,AxisString,Plot2D,CheckBox)
    XlblBefore = Plot1D.XLabel.String;
    XLim = Plot1D.XLim;
    YLim = Plot1D.YLim;
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
    
    if strcmpi(Plot1D.XLabel.String,XlblBefore) && CheckBox.Value
        Plot1D.XLim = XLim;
        Plot1D.YLim = YLim;
    end
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
