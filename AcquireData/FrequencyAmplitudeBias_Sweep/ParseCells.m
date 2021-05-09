function Output = ParseCells(Input, Parse)
if ~iscell(Input)
    UnCelledFlag = true;
    Input = {Input};
else
    UnCelledFlag = false;
end
Numel = numel(Input);

Output = cell(1,Numel);
j = 1;
for ic = 1:numel(Input)
    Output{j} = Parse(Input{ic});
    j = j+1;
end
if Numel==1 && UnCelledFlag
    Output = Output{:};
end
end

