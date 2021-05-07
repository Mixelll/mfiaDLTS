function Output = ParseStr2NumSym(Input)
if ~iscell(Input)
    Input = {Input};
end
Numel = 0;
for i = 1:numel(Input)
    if isstring(Input{i})
        Numel = Numel + numel(Input{i});
    else
        Numel = Numel + 1;
    end
end
    
StrChar = @(s) isstring(s) || ischar(s);
Output = cell(1,Numel);
j = 1;
for ic = Input
    in = ic{:};
    if StrChar(in) 
        if isstring(in) && numel(in)>1
            for s = in
                Output{j} = parse(s);
                j = j+1;
            end
        else
            Output{j} = parse(in);
            j = j+1;
        end
    else
        Output{j} = in;
        j = j+1;
    end
end
if Numel==1
    Output = Output{:};
end
function out = parse(in)
    if isempty(in)
        out = [];
    else
        Temp2 = str2num(in);
        if isempty(Temp2)
            out = sym(in);
        else
            out = Temp2;
        end
    end
end
end

