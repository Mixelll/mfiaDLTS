function out = parse_str2num_cell(in)
if isstring(in)
    in = char(in);
end
if isempty(in)
    out = [];
elseif in(1)=='{' && in(end)=='}'
    ins = in(2:end-1);
    if ~isempty(ins) && (any(regexp(ins,'([^\d\,\''](?=[^\'']))|((?<=[^\''])[^\d\,\''])')) || (~any(strfind(ins,',')) && ins(1)~='''' && ~any(str2num(ins))))
        out = strsplit(ins,',');
        for ic = 1:numel(out)
            Temp = str2num(out{ic});
            if ~isempty(Temp)
                out{ic} = Temp;
            end
        end  
    else
        out = eval(in);
    end
    ind = [];
    for ic = 1:numel(out)
        if isempty(out{ic})
            ind(end+1) = ic;
        end
    end
    out(ind) = [];
    if isempty(out), out = {}; end
else
    Temp2 = str2num(in);
    if isempty(Temp2)
        out = in;
    else
        out = Temp2;
    end
end
end