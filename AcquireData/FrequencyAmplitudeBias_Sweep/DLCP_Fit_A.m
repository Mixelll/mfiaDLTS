classdef DLCP_Fit_A < matlab.mixin.SetGet
	properties
        Range = [-inf inf]
        Area = (150e-4)^2 * pi
        RelativePermitivitty = 11.68
        Order = 1;
        FitPlotProperties = {};
        LegendProp = {'Location','Best'}
    end
	methods 
        function Fit(o, x,y, Target, varargin)
            [FitLeg, fun, span] = DLCP_fit(x,y,o.Area, 'range',o.Range, 'es',o.RelativePermitivitty, 'Order', o.Order);
            hold(Target, 'on')
            if ~isempty(varargin)
                legend(varargin{:});
            elseif isempty(isempty(Target.Legend))
                legend('Measured Data');
            end
            Lines = findobj(Target.Children, 'Type','Line');
            
            if ~isempty(Lines)
                fplot(fun, span, o.FitPlotProperties{:}, 'Parent',Target, 'Color',Lines(1).Color);
            else
                fplot(fun, span, o.FitPlotProperties{:}, 'Parent',Target);
            end
            if ~isempty(o.LegendProp)
                legend(o.LegendProp{:});
            end
            legend([Target.Legend.String(1:end-1) {[FitLeg newline 'Model: ' Target.Legend.String{end}]}])
            hold(Target, 'off')
        end
        function OK = Menu(o)
            OK = false;
            Properties = properties(o)';
            Values = get(o,Properties);
            answer = inputdlg(Properties,'Input Fit Parameters',1,ParseCells(Values,@parse_num_cell_sym2char));
            if ~isempty(answer)
                NewValues = ParseCells(answer,@parse_str2num_cell_sym);
                Pairs = [Properties ; NewValues];
                set(o, Pairs{:})
                OK = true;
            end
        end
	end
end