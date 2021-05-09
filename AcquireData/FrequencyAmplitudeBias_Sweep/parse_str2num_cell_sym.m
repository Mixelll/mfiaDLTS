function out = parse_str2num_cell_sym(in)
if isstring(in)
    in = char(in);
end
if isempty(in)
    out = [];
elseif in(1)=='{' && in(end)=='}'
%     out = strsplit(in(2:end-1),',');
    out = eval(in);
    ind = [];
    for ic = 1:numel(out)
        if isempty(out{ic})
            ind(end+1) = ic;
%             if out{ic}(1) == 
        end
    end
    out(ind) = [];
    if isempty(out), out = {}; end
else
    Temp2 = str2num(in);
    if isempty(Temp2)
        out = sym(in);
    else
        out = Temp2;
    end
end
end