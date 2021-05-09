function Output = ParseCellsAndStringArrays(Input, Parse)
if ~iscell(Input)
    UnCelledFlag = true;
    Input = {Input};
else
    UnCelledFlag = false;
end
Numel = 0;
for i = 1:numel(Input)
    if isstring(Input{i})
        Numel = Numel + numel(Input{i});
    else
        Numel = Numel + 1;
    end
end

Output = cell(1,Numel);
j = 1;
for ic = 1:numel(Input)
    in = Input{ic};
    if isstring(in) && numel(in)>1
        for s = in
            Output{j} = Parse(s);
            j = j+1;
        end
    else
        Output{j} = Parse(in);
        j = j+1;
    end
end
if Numel==1 && UnCelledFlag
    Output = Output{:};
end
end

