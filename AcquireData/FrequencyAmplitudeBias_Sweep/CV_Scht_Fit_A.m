classdef CV_Scht_Fit_A < matlab.mixin.SetGet
	properties
        Range = [-inf inf]
        Area = (150e-4)^2 * pi
        RelativePermitivitty = 11.68
        RelativePermitivitty_Limits = [-inf 0 inf]
        Doping = sym('N')
        Doping_Limits = [-inf 0 inf]
        Vb = sym('Vb')
        Vb_Limits = [-inf 0 inf]
        IdealityFactor = sym('n')
        IdealityFactor_Limits = [-inf 0 inf]
        FitProperties = {}
        FitPlotProperties = {}
    end
	methods 
        function Fit(o, x,y, Target, varargin)
            LimChanged = @(x) any(x~=[-inf 0 inf]);
            Limits = {};
            if LimChanged(o.RelativePermitivitty_Limits), Limits(end+1:end+2) = {'es' o.RelativePermitivitty_Limits}; end
            if LimChanged(o.Doping_Limits), Limits(end+1:end+2) = {'N' o.Doping_Limits}; end
            if LimChanged(o.Vb_Limits), Limits(end+1:end+2) = {'Vb' o.Vb_Limits}; end
            if LimChanged(o.IdealityFactor_Limits), Limits(end+1:end+2) = {'n' o.IdealityFactor_Limits}; end
            [FitLeg, fun, span] = C_schot_fit_A(x,y,o.Range,o.Area,o.RelativePermitivitty,o.Doping,o.Vb,o.IdealityFactor,o.FitProperties,Limits{:});
            hold(Target, 'on')
            if ~isempty(varargin)
                if iscell(varargin{:})
                   legend([varargin{:}(1:end-1) {[varargin{:}{end} newline FitLeg]}]);
                else
                    legend([char(varargin{:}) newline FitLeg]);
                end
            else
                Leg0 = legend(FitLeg);
            end
            Lines = findobj(Target.Children, 'Type','Line');
            
            if ~isempty(Lines)
                fplot(fun, span, o.FitPlotProperties{:}, 'Parent',Target, 'Color',Lines(1).Color);
            else
                fplot(fun, span, o.FitPlotProperties{:}, 'Parent',Target);
            end
            hold(Target, 'off')
        end
        function Menu(o)
            Properties = properties(o)';
            Values = get(o,Properties);
            answer = inputdlg(Properties,'Input Fit Parameters',1,ParseCells(Values,@parse_num_cell_sym2char));
            if ~isempty(answer)
                NewValues = ParseCells(answer,@parse_str2num_cell_sym);
                Pairs = [Properties ; NewValues];
                set(o, Pairs{:})
            end
        end
	end
end