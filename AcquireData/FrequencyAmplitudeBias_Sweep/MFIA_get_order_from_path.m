function order = MFIA_get_order_from_path(path, varargin)
if isempty(varargin)
    search_order = 'order_'; % input the string just before the first parameter
else
    search_order = varargin{:};
end
order_string = path((strfind(path,search_order)+length(search_order)):end);
order_ind = strfind(order_string, '_');
order = {order_string(1:order_ind(1)-1), order_string(order_ind(1)+1:order_ind(2)-1), order_string(order_ind(2)+1:end)};
end