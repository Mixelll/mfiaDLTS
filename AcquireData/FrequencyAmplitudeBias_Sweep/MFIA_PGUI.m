function [figs] = MFIA_PGUI(SubPlots,D)
clear 'DataFigureHandles'
k = 0;
Nsbp = numel(SubPlots);
if D==3
    for j = 1:size(SubPlots,2)
        for i = 1:size(SubPlots,1)
            if strcmpi(class(SubPlots(i,j)),'matlab.graphics.axis.Axes')
            k = k+1;
%             FigUserData.type = SubPlots(i,j).Title.UserData;
%             fig = figure('Position',get(0,'Screensize'), 'UserData',FigUserData);
            fig = figure('Position',get(0,'Screensize'));
            AddProperties(fig);
            ax3D = subplot(1,3,1);
            ax2D = subplot(1,3,2); 
            ax1D = subplot(1,3,3);

            copyobj(get(SubPlots(i,j), 'Children'),ax3D)
            ax3D = update_structure(ax3D, SubPlots(i,j), 'ignore',{'Parent'}, 'new', true);
            fig.UserData.type = ax3D.UserData.type;
            ax2D.UserData.type = ax3D.UserData.type;
            ax1D.UserData.type = ax3D.UserData.type;
            if Nsbp>2 || Nsbp==1
                ax3D.Position = [0.05    0.11    0.304    0.815];
            else
                ax3D.Position = [0.05 ax3D.Position(2:end)];
            end
            colorbar(ax3D, 'north')
            
            P3D = ax3D.Position; P2D = ax2D.Position; P1D = ax1D.Position;
            Temp = ax3D.UserData.axesvec;
            Tempf = fieldnames(Temp)';
            ax_lbls = {Temp.(Tempf{1}).label Temp.(Tempf{2}).label Temp.(Tempf{3}).label};
            ax_ranges = {[min(Temp.(Tempf{1}).vec) max(Temp.(Tempf{1}).vec)] [min(Temp.(Tempf{2}).vec) max(Temp.(Tempf{2}).vec)] [min(Temp.(Tempf{3}).vec) max(Temp.(Tempf{3}).vec)]};
            ax_lengths = {length(Temp.(Tempf{1}).vec) length(Temp.(Tempf{2}).vec) length(Temp.(Tempf{3}).vec)};
            PopStrVal = [ax_lbls ; ax_ranges ; Tempf; ax_lengths];
            PopOne2D = PopStrVal{2,2};
            PopOne3D = PopStrVal{2,1};
            SliderStep2 = max(1/PopStrVal{4,2},0.1);

            
            c1D_Lim = uicontrol('Parent',fig,'Style','checkbox', 'Units','normalized', 'value',0, 'max',1, 'min',0, 'String','Keep Axes Lim', 'FontSize',13, 'Position',[P1D(1) P1D(2)-0.07 0.08 0.02]);
            c1D_Hold = uicontrol('Parent',fig,'Style','checkbox', 'Units','normalized', 'value',0, 'max',1, 'min',0, 'String','Hold Plot', 'FontSize',13, 'Position',[P1D(1)+0.1 P1D(2)-0.07 0.08 0.02]);
            c2D_Lim = uicontrol('Parent',fig,'Style','checkbox', 'Units','normalized', 'value',0, 'max',1, 'min',0, 'String','Keep Axes Lim', 'FontSize',13, 'Position',[P2D(1) P2D(2)-0.07 0.08 0.02]);
            
            PopCall1D = @(src,cbdata) popmenu_callbackFit1D(src,cbdata, ax1D);
            pf1D = uicontrol('Parent',fig, 'Style','popupmenu', 'string',[{'Disable'} RegisteredNames('Fit Classes')], 'FontSize',10, 'Units','normalized', 'Position', [P1D(1) P1D(2)-0.1 0.1 0.03], 'Callback',PopCall1D);
            
            SliderCallFromPop2D = @(src,cbdata,AxisString) slider_callback2D(src,cbdata,ax1D,AxisString,ax2D,c1D_Lim,c1D_Hold,pf1D);
            s2D = uicontrol('Parent',fig,'Style','slider', 'Units','normalized', 'Callback', @(src,cbdta) SliderCallFromPop2D(src,cbdta,PopStrVal{3,2}),...
                'value',mean(PopOne2D), 'min',PopOne2D(1), 'max',PopOne2D(2), 'SliderStep',[1/PopStrVal{4,2} SliderStep2]);
            s2D.Position = [P2D(1) P2D(2)-0.11 P2D(3) 0.025];
            LeftTxt2D = uicontrol('Style','text', 'Units','normalized', 'Position',[P2D(1)-0.03 P2D(2)-0.11 0.03 0.025],'String',num2str(PopOne2D(1),3) ,'FontSize',15);
            RightTxt2D = uicontrol('Style','text', 'Units','normalized', 'Position',[P2D(1)+P2D(3) P2D(2)-0.11 0.03 0.025],'String',num2str(PopOne2D(2),3) ,'FontSize',15);

            PopCall2D = @(src,cbdata,PopStrVal) popmenu_callback2D(src,cbdata,s2D,SliderCallFromPop2D,PopStrVal,LeftTxt2D,RightTxt2D,pf1D);
            p2D = uicontrol('Parent',fig, 'Style','popupmenu', 'string',PopStrVal(1,2:3), 'FontSize',10, 'Units','normalized', 'Position', [P2D(1)+0.1 P2D(2)-0.088 0.13 0.03], 'Callback',@(src,cbdta)PopCall2D(src,cbdta,PopStrVal(:,2:3)));

            SliderCallFromPop3D = @(src,cbdata,AxisString) slider_callback3D(src,cbdata,ax2D,AxisString,ax3D.UserData,c2D_Lim);
            s3D = uicontrol('Parent',fig,'Style','slider', 'Units','normalized', 'Callback', @(src,cbdta) SliderCallFromPop3D(src,cbdta,PopStrVal{3,1}),...
                'value',mean(PopOne3D), 'min',PopOne3D(1), 'max',PopOne3D(2), 'SliderStep',[1/PopStrVal{4,1} SliderStep2]);
            s3D.Position = [P3D(1) P3D(2)-0.11 P3D(3) 0.025];
            LeftTxt3D = uicontrol('Style','text', 'Units','normalized', 'Position',[P3D(1)-0.03 P3D(2)-0.11 0.03 0.025],'String',num2str(PopOne3D(1),3) ,'FontSize',15);
            RightTxt3D = uicontrol('Style','text', 'Units','normalized', 'Position',[P3D(1)+P3D(3) P3D(2)-0.11 0.03 0.025],'String',num2str(PopOne3D(2),3) ,'FontSize',15);

            PopCall3D = @(src,cbdata) popmenu_callback3D(src,cbdata,s3D,SliderCallFromPop3D,p2D,PopCall2D,PopStrVal,LeftTxt3D,RightTxt3D);
            %p3D =
            uicontrol('Parent',fig, 'Style','popupmenu', 'string',PopStrVal(1,:), 'FontSize',10, 'Units','normalized', 'Position', [P3D(1)+0.13 P3D(2)-0.088 0.13 0.03], 'Callback',PopCall3D);
            
            
            ExportPlotBtnCall = @(src,cbdata) ExportPlotCallBack(src, cbdata, ax1D);
            %ExportPlotBtn = 
            uicontrol('Parent',fig, 'string','Export Plot', 'FontSize',10, 'Units','normalized', 'Position', [P1D(1)+0.105 P1D(2)-0.1 0.05 0.03], 'Callback',ExportPlotBtnCall);

            %FigPropBtn = 
            uicontrol('Parent',fig, 'string','Settings', 'FontSize',10, 'Units','normalized', 'Position', [P1D(1)+0.16 P1D(2)-0.1 0.05 0.03], 'Callback',@ChangeFigProperties);
            
            TextCall1D = @(src,cbdata) textbox_callback2D(src,cbdata, s2D);
            %SliceText =
            uicontrol('Parent',fig, 'Style','edit', 'string','0', 'FontSize',10, 'Units','normalized', 'Position', [P2D(1)+0.06 P2D(2)-0.08 0.03 0.03], 'Callback',TextCall1D);
            
            figs(k) = fig;
            end
        end
    end
elseif D==2
    for j = 1:size(SubPlots,2)
        for i = 1:size(SubPlots,1)
            if strcmpi(class(SubPlots(i,j)),'matlab.graphics.axis.Axes')
            k = k+1;
            fig = figure('Position',get(0,'Screensize'));
            AddProperties(fig);
            ax2D = subplot(1,3,1);
            ax2Dpc = subplot(1,3,2); 
            ax1D = subplot(1,3,3);

            copyobj(get(SubPlots(i,j), 'Children'),ax2D)
            ax2D = update_structure(ax2D, SubPlots(i,j), 'ignore',{'Parent'}, 'new', true);
            fig.UserData.type = ax2D.UserData.type;
            ax2Dpc.UserData.type = ax2D.UserData.type;
            ax1D.UserData.type = ax2D.UserData.type;
            LimCall = @(src, evt) LimCallback2D(src, evt, ax2D, ax2Dpc);
            ax2D.CreateFcn = @(src, evt) addlistener(ax2D, 'XLim', 'PostSet', LimCall);
            addlistener(ax2D, 'XLim', 'PostSet', LimCall);
 
            try 
                pcolor(ax2Dpc, mean(ax2D.Children.XData,2), mean(ax2D.Children.YData,1), ax2D.Children.CData.');
            catch
                pcolor(ax2Dpc, mean(ax2D.Children.XData,2), mean(ax2D.Children.YData,1), ax2D.Children.CData);
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
            P2D = ax2D.Position; P2Dpc = ax2Dpc.Position; P1D = ax1D.Position;
            Temp = ax2D.UserData.axesvec;
            Tempf = fieldnames(Temp)';
            for ic = 1:length(Tempf)
                ax_lbls_all{ic} = Temp.(Tempf{ic}).label;
            end
            Tempf = [Tempf(strcmpi(ax_lbls_all, ax2D.XLabel.String)) Tempf(strcmpi(ax_lbls_all, ax2D.YLabel.String))];
            ax_lbls = {Temp.(Tempf{1}).label Temp.(Tempf{2}).label};
            ax_ranges = {[min(Temp.(Tempf{1}).vec) max(Temp.(Tempf{1}).vec)] [min(Temp.(Tempf{2}).vec) max(Temp.(Tempf{2}).vec)]};
            ax_lengths = {length(Temp.(Tempf{1}).vec) length(Temp.(Tempf{2}).vec)};
            PopStrVal = [ax_lbls ; ax_ranges ; Tempf; ax_lengths];
            PopOne2D = PopStrVal{2,1};
            SliderStep2 = max(1/PopStrVal{4,2},0.1);
            
            c1D_Lim = uicontrol('Parent',fig,'Style','checkbox', 'Units','normalized', 'value',0, 'max',1, 'min',0, 'String','Keep Axes Lim', 'FontSize',13, 'Position',[P1D(1) P1D(2)-0.07 0.08 0.02]);
            c1D_Hold = uicontrol('Parent',fig,'Style','checkbox', 'Units','normalized', 'value',0, 'max',1, 'min',0, 'String','Hold Plot', 'FontSize',13, 'Position',[P1D(1)+0.1 P1D(2)-0.07 0.08 0.02]);
            
            PopCall1D = @(src,cbdata) popmenu_callbackFit1D(src,cbdata, ax1D);
            pf1D = uicontrol('Parent',fig, 'Style','popupmenu', 'string',[{'Disable'} RegisteredNames('Fit Classes')], 'FontSize',10, 'Units','normalized', 'Position', [P1D(1) P1D(2)-0.1 0.1 0.03], 'Callback',PopCall1D);
            
            SliderCallFromPop2D = @(src,cbdata,AxisString) slider_callback2D(src,cbdata,ax1D,AxisString,ax2Dpc,c1D_Lim,c1D_Hold,pf1D);
            s2D = uicontrol('Parent',fig,'Style','slider', 'Units','normalized', 'Callback', @(src,cbdta) SliderCallFromPop2D(src,cbdta,PopStrVal{3,2}),...
                'value',mean(PopOne2D), 'min',PopOne2D(1), 'max',PopOne2D(2), 'SliderStep',[1/PopStrVal{4,2} SliderStep2]);
            s2D.Position = [P2Dpc(1) P2Dpc(2)-0.11 P2Dpc(3) 0.025];
            LeftTxt2D = uicontrol('Style','text', 'Units','normalized', 'Position',[P2Dpc(1)-0.03 P2Dpc(2)-0.11 0.03 0.025],'String',num2str(PopOne2D(1),3) ,'FontSize',15);
            RightTxt2D = uicontrol('Style','text', 'Units','normalized', 'Position',[P2Dpc(1)+P2Dpc(3) P2Dpc(2)-0.11 0.03 0.025],'String',num2str(PopOne2D(2),3) ,'FontSize',15);

            PopCall2D = @(src,cbdata,PopStrVal) popmenu_callback2D(src,cbdata,s2D,SliderCallFromPop2D,PopStrVal,LeftTxt2D,RightTxt2D,pf1D);
            p2D = uicontrol('Parent',fig, 'Style','popupmenu', 'string',PopStrVal(1,:), 'FontSize',10, 'Units','normalized', 'Position', [P2Dpc(1)+0.1 P2Dpc(2)-0.088 0.13 0.03], 'Callback',@(src,cbdta)PopCall2D(src,cbdta,PopStrVal(:,:)));

            ExportPlotBtnCall = @(src,cbdata) ExportPlotCallBack(src, cbdata, ax1D);
            %ExportPlotBtn = 
            uicontrol('Parent',fig, 'string','Export Plot', 'FontSize',10, 'Units','normalized', 'Position', [P1D(1)+0.105 P1D(2)-0.1 0.05 0.03], 'Callback',ExportPlotBtnCall);
            
            TextCall1D = @(src,cbdata) textbox_callback2D(src,evt, s2D);
            %FigPropBtn = 
            uicontrol('Parent',fig, 'string','Settings', 'FontSize',10, 'Units','normalized', 'Position', [P1D(1)+0.16 P1D(2)-0.1 0.05 0.03], 'Callback',@ChangeFigProperties);
            
            TextCall1D = @(src,cbdata) textbox_callback2D(src,cbdata, s2D);
            %SliceText =
            uicontrol('Parent',fig, 'Style','edit', 'string','0', 'FontSize',10, 'Units','normalized', 'Position', [P2Dpc(1)+0.06 P2Dpc(2)-0.08 0.03 0.03], 'Callback',TextCall1D);
            
            figs(k) = fig;
            end
        end
    end
end
%% Functions
function slider_callback3D(src,~,Plot2D,SelectAxString,DataStruct,CheckBoxHold)
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
        
    pcolor(Plot2D, Axes.(AxNames{2}).vec, Axes.(AxNames{1}).vec, C);
    colormap(Plot2D, 'turbo')
    colorbar(Plot2D, 'northoutside')
    Plot2D.XLabel.String = Axes.(AxNames{2}).label;
    Plot2D.YLabel.String = Axes.(AxNames{1}).label;
    Plot2D.Title.String =[Axes.(SelectAxString).label ' = ' num2str(Axes.(SelectAxString).vec(SelectMinInd))];
    Plot2D.UserData = UserData;
    if strcmpi(Plot2D.XLabel.String,XlblBefore) && strcmpi(Plot2D.YLabel.String,YlblBefore) && CheckBoxHold.Value
        Plot2D.XLim = XLim;
        Plot2D.YLim = YLim;
    end  
end

function slider_callback2D(src,~,Plot1D,AxisString,Plot2D,CheckBoxLim,CheckBoxHold,pf1D)
    Hold = CheckBoxHold.Value;
    TitleBefore = Plot1D.Title.String;
    if ~isempty(Plot1D.Legend)
        LegendBefore = Plot1D.Legend.String;
    else
        LegendBefore = {};
    end
    XlblBefore = Plot1D.XLabel.String;
    XLim = Plot1D.XLim;
    YLim = Plot1D.YLim;
    if ~CheckBoxLim.Value
        Plot1D.XAxis.LimitsMode = 'auto';
        Plot1D.YAxis.LimitsMode = 'auto';
    end
    UserData = Plot1D.UserData;
    P2DC = Plot2D.Children;
    P2DX = Plot2D.XLabel.String;
    P2DY = Plot2D.YLabel.String;
    PlotFlag = 0;
    if strcmpi(P2DX,AxisString)
        [~,SelectMinInd] = min(abs(P2DC.XData - src.Value));
        x = P2DC.YData;
        y = shiftdim(P2DC.CData(:,SelectMinInd));
        Title = [P2DX ' = ' num2str(P2DC.XData(SelectMinInd))];
        XLabel = P2DY;
        PlotFlag = 1;
    elseif strcmpi(P2DY,AxisString)
        [~,SelectMinInd] = min(abs(P2DC.YData - src.Value));
        x = P2DC.XData;
        y = shiftdim(P2DC.CData(SelectMinInd,:));
        Title = [P2DY ' = ' num2str(P2DC.YData(SelectMinInd))];
        XLabel = P2DX;
        PlotFlag = 1;
    end
    
    if Hold && strcmpi(XlblBefore, XLabel)
        hold(Plot1D, 'on')
    else
        Hold = false;
    end
    if PlotFlag
        plot(Plot1D, x, y, Plot1D.Parent.UserData.Properties.Plot_Style)
        Plot1D.Title.String = Title;
        Plot1D.XLabel.String = XLabel;
        Plot1D.UserData = UserData;
        Plot1D.YLabel.String = Plot1D.UserData.type;
        if Hold && (numel(findobj(Plot1D.Children, 'Type','Line'))>2 || ~isempty(LegendBefore))
            legend(Plot1D, [LegendBefore Plot1D.Title.String])
        elseif Hold
            legend(Plot1D, {TitleBefore Plot1D.Title.String});
        end
        LegFitIn = Title;
        if ~isempty(Plot1D.Legend)
            LegFitIn = Plot1D.Legend.String;
        end
        if strcmpi(Plot1D.XLabel.String,XlblBefore) && CheckBoxLim.Value
            Plot1D.XLim = XLim;
            Plot1D.YLim = YLim;
        end
        if ~strcmpi(pf1D.String(pf1D.Value), 'Disable')
            Plot1D.Parent.UserData.Fit.(pf1D.String{pf1D.Value}).Fit(x,y, Plot1D, LegFitIn)
        end
        hold(Plot1D, 'off')
    end
end
function textbox_callback2D(src,~, slider2D)
    srcc.Value = str2double(src.String);
    if ~isnan(srcc.Value)
        slider2D.Callback(srcc,srcc);
    end
end
function popmenu_callback3D(src,~,slider3D,SliderCall,Pop2D,Pop2D_call,value_cell,LeftTxt,RightTxt)
    SV = src.Value;
    AxBounds = value_cell{2,SV};
    set(slider3D, 'value',mean(AxBounds), 'min',AxBounds(1), 'max',AxBounds(2), 'SliderStep',[1/value_cell{4,SV} max(1/value_cell{4,SV},0.1)],...
        'Callback', @(src,cbdata) SliderCall(src,cbdata,value_cell{3,SV}))
    set(LeftTxt, 'String', num2str(AxBounds(1),3))
    set(RightTxt, 'String', num2str(AxBounds(2),3))
    value_cell(:,SV) = [];
    set(Pop2D, 'string',value_cell(1,:), 'Callback',@(src,cbdata)Pop2D_call(src,cbdata,value_cell));
end

function popmenu_callback2D(src,~,slider2D,SliderCall,value_cell,LeftTxt,RightTxt,pf1D)
    SV = src.Value;
    AxBounds = value_cell{2,SV};
    set(slider2D, 'value',mean(AxBounds), 'min',AxBounds(1), 'max',AxBounds(2),  'SliderStep',[1/value_cell{4,SV} max(1/value_cell{4,SV},0.1)],...
        'Callback', @(src,cbdata) SliderCall(src,cbdata,value_cell{1,SV}))
    set(LeftTxt, 'String', num2str(AxBounds(1),3))
    set(RightTxt, 'String', num2str(AxBounds(2),3))
    set(pf1D, 'value', 1)
end

function popmenu_callbackFit1D(src,~, Plot1D)
    SV = src.Value;
    pop = src;
    SF = pop.String{SV};
    [FitNames,FitHandles] = RegisteredNames('Fit Classes');
    PopFitNames = [{'Disable'} FitNames];
    if length(unique(lower(PopFitNames)))~=length(pop.String) || ~all(strcmpi(unique(lower(PopFitNames)), pop.String'))
        pop.String = PopFitNames;
    end
    if ~strcmpi(SF, 'Disable')
        FC = FitHandles{strcmpi(FitNames, SF)};
        ParentFig = pop.Parent;
        if isfield(ParentFig.UserData, 'Fit')
            ExFitNames = fieldnames(ParentFig.UserData.Fit);
        else
            ExFitNames = {};
        end
        if ~any(strcmp(ExFitNames, SF))
            ParentFig.UserData.Fit.(SF) = FC();
        end
        OK = ParentFig.UserData.Fit.(SF).Menu();
        Lines = findobj(Plot1D.Children, 'Type','Line');
        Line = Lines(1);
        LegFitIn = Plot1D.Title.String;
        if ~isempty(Plot1D.Legend)
            LegFitIn = Plot1D.Legend.String;
        end
        if OK
            ParentFig.UserData.Fit.(SF).Fit(Line.XData,Line.YData, Plot1D, LegFitIn)
        end
    end
end

function LimCallback2D(~, ~, Plot2Dsurf, Plot2Dpc)
    XLim = Plot2Dsurf.XLim;
    YLim = Plot2Dsurf.YLim;
    ZLim = Plot2Dsurf.ZLim;
    XVec = mean(Plot2Dsurf.Children.XData,2);
    YVec = mean(Plot2Dsurf.Children.YData,1);
    CData = Plot2Dsurf.Children.CData;

    CData(CData<ZLim(1) | CData>ZLim(2)) = NaN;
    XVecL = XVec>XLim(1) & XVec<XLim(2);
    YVecL = YVec>YLim(1) & YVec<YLim(2);
    XVec = XVec(XVecL);
    YVec = YVec(YVecL);
    CData = CData(XVecL,YVecL);
    try 
        pcolor(Plot2Dpc, XVec, YVec, CData.');
    catch
        pcolor(Plot2Dpc, XVec, YVec, CData);
    end
    colorbar(Plot2Dpc, 'northoutside')
    Plot2Dpc.XLabel.String = Plot2Dsurf.XLabel.String;
    Plot2Dpc.YLabel.String = Plot2Dsurf.YLabel.String;
    Plot2Dpc.UserData.type = Plot2Dsurf.UserData.type;
end

function ExportPlotCallBack(~, ~, Plot1D)
    prompt = {'Enter plot number:', 'Enter legend:', 'Enter font size:', 'Enter plot properies:', 'Enter Range:', 'Copy As is:'};
    dlgtitle = 'Add Plot to Target Figure';
    [Plot1D.Parent.UserData.ExportPlotProperties, OK] = StructrureFieldsMenu(Plot1D.Parent.UserData.ExportPlotProperties, @parse_num_cell_sym2char, @parse_str2num_cell, prompt, dlgtitle);
    EPP = Plot1D.Parent.UserData.ExportPlotProperties;
    if OK && ~isempty(EPP.Plot_Number)
        Legend = EPP.Legend;
        if isempty(Legend)
            if ~isempty(Plot1D.Legend)
                Legend = Plot1D.Legend.String;
            elseif ~isempty(Plot1D.Title)
                Legend = Plot1D.Title.String;
            end
        end
        Range = EPP.Range;
        if length(Range)<2
            Lx = -inf; Hx = -inf;
        else
            Lx = Range(1); Hx = Range(2);
        end
        CopyChildren = EPP.CopyChildren;
        
        xlbl = Plot1D.XLabel.String;
        ylbl = Plot1D.YLabel.String;

        Children = Plot1D.Children(end:-1:1);
        Cells = cell(1,numel(Children));
        for i2 = 1:numel(Children)
            Children(i2).UserData = [];
            xd = Children(i2).XData;
            if strcmpi(Children(i2).Type, 'line') && ~CopyChildren
                Cells{i2} = [Children(i2).XData(xd>=Lx & xd<=Hx); Children(i2).YData(xd>=Lx & xd<=Hx)];
            else
                if strcmpi(Children(i2).Type, 'functionline') && ~CopyChildren
                    xr = Children(i2).XRange;
                    Children(i2).UserData = ['XRange=[' num2str(max(Lx,xr(1))) ' ' num2str(min(Hx,xr(2))) ']'];
                end
                Cells{i2} = Children(i2);
            end
        end
        s_plot([], Cells, EPP.PlotProperties, '', '', Legend, EPP.Plot_Number, xlbl, ylbl, EPP.FontSize, '', 0, 0, 1);
    end
end

function ChangeFigProperties(src, ~)
    AddProperties(src.Parent);
    src.Parent.UserData.Properties = StructrureFieldsMenu(src.Parent.UserData.Properties, @parse_num_cell_sym2char, @parse_str2num_cell);
end
function AddProperties(fig)
    % General properties
    MoreProp.Plot_Style = 'o';
    if isfield(fig.UserData, 'Properties')
        fig.UserData.Properties = update_structure(fig.UserData.Properties, MoreProp, 'onlynew',1);
    else
        fig.UserData.Properties = MoreProp;
    end
    % Export plot properties
    MoreExportPlotProp.Plot_Number = '';
    MoreExportPlotProp.Legend = {};
    MoreExportPlotProp.FontSize = '';
    MoreExportPlotProp.PlotProperties = {};
    MoreExportPlotProp.Range = [-inf inf];
    MoreExportPlotProp.CopyChildren = 1;
    if isfield(fig.UserData, 'AddPlotProperties')
        fig.UserData.ExportPlotProperties = update_structure(fig.UserData.ExportPlotProperties, MoreExportPlotProp, 'onlynew',1);
    else
        fig.UserData.ExportPlotProperties = MoreExportPlotProp;
    end
end
end

