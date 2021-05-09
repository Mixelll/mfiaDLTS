function Output = ParseNumSym2Str(Input)
NumSym = @(x) isnumeric(x) || isa(x,'sym');
StrChar = @(s) isstring(s) || ischar(s);
if ~iscell(Input)
    UnCelledFlag = true;
    Input = {Input};
else
    UnCelledFlag = false;
end
Numel = numel(Input);

Output = cell(1,Numel);
j = 1;
for ic = Input
    in = ic{:};
    if NumSym(in)
        Output{j} = parse_num_or_sym2str(in);
        j = j+1;
    else
        Output{j} = in;
        j = j+1;
    end
end
if Numel==1 && UnCelledFlag
    Output = Output{:};
end
%%
function out = parse_str2num_or_sym(in)
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
function out = parse_num_or_sym2str(in)
    if isempty(in)
        out = '';
    elseif isnumeric(in)
        out = num2str(in);
    else
        out = char(in);
    end
end
end

