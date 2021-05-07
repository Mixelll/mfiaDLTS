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
        function Fit(o, x,y, Target)
            LimChanged = @(x) any(x~=[-inf 0 inf]);
            Limits = {};
            if LimChanged(o.RelativePermitivitty_Limits), Limits(end+1:end+2) = {'es' o.RelativePermitivitty_Limits}; end
            if LimChanged(o.Doping_Limits), Limits(end+1:end+2) = {'N' o.Doping_Limits}; end
            if LimChanged(o.Vb_Limits), Limits(end+1:end+2) = {'Vb' o.Vb_Limits}; end
            if LimChanged(o.IdealityFactor_Limits), Limits(end+1:end+2) = {'n' o.IdealityFactor_Limits}; end
            [~, fun, span] = C_schot_fit_A(x,y,o.Range,o.Area,o.RelativePermitivitty,o.Doping,o.Vb,o.IdealityFactor,o.FitProperties,Limits{:});
            hold(Target, 'on')
            fplot(fun, span, o.FitPlotProperties{:}, 'Parent',Target);
            hold(Target, 'off')
        end
        function Menu(o)
            Properties = properties(o)';
            Values = get(o,Properties);
            Values
            ParseStr2NumSym(Values)
%             answer = inputdlg(Properties,'Input Fit Parameters',1,ParseStr2NumSym(Values))
        end
	end
end