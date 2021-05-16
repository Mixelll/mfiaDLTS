function [StrCells, HandleCells]= RegisteredNames(varargin)
StrCells = {};
HandleCells = {};
for c = varargin
    switch c{:}
    case 'Fit Classes'
        Str = {'CV_Scht_Fit_A', 'DLCP_Fit_A'};
        Handle = {@CV_Scht_Fit_A, @DLCP_Fit_A};
    end
    StrCells = [StrCells Str];
    HandleCells = [HandleCells Handle];
end
[StrCells, Indices] = unique(StrCells);
HandleCells = HandleCells(Indices);
end

