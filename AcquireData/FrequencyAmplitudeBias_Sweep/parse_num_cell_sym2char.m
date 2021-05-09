function out = parse_num_cell_sym2char(in)
if isstring(in)
    in = char(in);
end
if isempty(in)
    out = '';
elseif isnumeric(in)
    out = num2str(in);
elseif in(1)=='{' && in(end)=='}'
    
else
    try
        out = char(in);
    catch
        out = in;
    end
end
end