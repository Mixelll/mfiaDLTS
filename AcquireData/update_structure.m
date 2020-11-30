function [original_struct] = update_structure_exact(original_struct,updated_struct)
original_cell = fn_struct2cell(original_struct);
updated_cell = fn_struct2cell(updated_struct);
for c = original_cell
    field_cmp = strcmp(updated_cell(3,:), c{3});
    if any(field_cmp)
    eval(['original_struct' c{3} '=' updated_cell{1,field_cmp} ';'])
    end
end   
end

