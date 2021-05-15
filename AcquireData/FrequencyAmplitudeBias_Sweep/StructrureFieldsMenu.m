function OUT = StructrureFieldsMenu(p)
OUT = struct;
Cells = fn_struct2cell(p);
answer = inputdlg(Cells(1,:),'Input Parameters',1,ParseCells(Cells(2,:),@parse_num_cell_sym2char));
if ~isempty(answer)
    NewValues = ParseCells(answer,@parse_str2num_cell);
    Cells(2,:) = NewValues;
    for c = Cells
        Temp = c{2};
        eval(['OUT' c{3} '=Temp;'])
    end
end
end

