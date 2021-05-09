function out = parse_str2num_sym(in)
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