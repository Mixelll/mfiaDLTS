function out = parse_num_sym2char(in)
    if isempty(in)
        out = '';
    elseif isnumeric(in)
        out = num2str(in);
    else
        try
            out = char(in);
        catch
            out = in;
        end
    end
end