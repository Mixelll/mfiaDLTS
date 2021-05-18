function [OUT, OK] = StructrureFieldsMenu(p, ParseInMenu, ParseOutMenu, varargin)
OK = true;
OUT = struct;
Cells = fn_struct2cell(p);
Title = 'Input Parameters';
Prompt = Cells(1,:);
for c = varargin
    if ischar(c{:}) || isstring(c{:})
        Title = c{:};
    elseif  iscell(c{:})
        Prompt = c{:};
    end
end
answer = inputdlg(Prompt,Title,1,ParseCells(Cells(2,:),ParseInMenu));
if ~isempty(answer)
    NewValues = ParseCells(answer,ParseOutMenu);
    Cells(2,:) = NewValues;
    for c = Cells
        Temp = c{2};
        eval(['OUT' c{3} '=Temp;'])
    end
else
    OK = false;
end
end

