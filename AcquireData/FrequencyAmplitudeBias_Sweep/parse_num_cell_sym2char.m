function out = parse_num_cell_sym2char(in)
StrChar = @(s) isstring(s) || ischar(s);
if isstring(in)
    in = char(in);
end
if isempty(in) && ~iscell(in)
    out = '';
elseif isnumeric(in)
    out = num2str(in);
% elseif in(1)=='{' && in(end)=='}'
elseif iscell(in)
    out = '{';
    for ic = 1:numel(in)
        if StrChar(in{ic})
            temp = ['''' char(in{ic}) ''''];
        else
            temp = parse_num_sym2char(in{ic});
        end
            out = [out temp ','];
    end
    if out(end)==',', out(end)='}'; else, out(end+1)='}'; end
else
    try
        out = char(in);
    catch
        out = in;
    end
end
end