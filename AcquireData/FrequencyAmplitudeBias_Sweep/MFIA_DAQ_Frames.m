function [DataTable, t_delta] = MFIA_DAQ_Frames(DataCells, varargin)
CellStrChar = @(s) iscell(s) || isstring(s) || ischar(s);
ChrStr = @(s) ischar(s) || isstring(s);
CellStrCharFn = @(s) iscell(s) || isstring(s) || ischar(s) || isa(s, 'function_handle');
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('DataSelect', {'Capacitance'}, CellStrChar); % Select data type
p.addParameter('TemperatureRange', [-inf inf 1], @isnumeric); % Select temperature range. Third element: true to select inside range, false for outside 
p.addParameter('SetNumbers', [], @isnumeric); % Select sets by numbers
p.addParameter('Range', [-inf inf], @isnumeric); % Trim data by time range

p.parse(varargin{:});

DS = p.Results.DataSelect;
TR = p.Results.TemperatureRange;
SN = p.Results.SetNumbers;
R = p.Results.Range;
if isempty(R)
    R = [-inf inf -inf inf];
elseif length(R)<4
    R = [R(1:2) -inf inf];
end
if isempty(TR)
    TR = [-inf inf 1];
end
if length(TR)<3
    TeRflag = true;
else
    TeRflag = logical(TR(3));
end

Ncells = size(DataCells,2);
Fnames = fieldnames(DataCells{2,1})';
DataTable = table('Size',[0,9],'VariableNames', {'Data Type', 'Set Number', 'Set Name', 'T', 'Time Constant', 't0', 'Length Seconds', 'Length', 'Data'},...
    'VariableTypes', {'string','double','string','double','double','double','double','double','cell'});
x_delta_vec = zeros(1,Ncells);
for i = 1:Ncells
    Data = DataCells{2,i};
    xV = Data.(Fnames{1})(:,1); 
    x_delta = mean(diff(xV));
    
    for c = Fnames
        DataC = Data.(c{:});
        DataTable(end+1,:) = [c {i} DataCells(1,i) DataCells{3,i}(2) x_delta DataC(1,1) DataC(end,1) length(DataC) {DataC}];
    end
    x_delta_vec(i) = x_delta;
end
t_delta = mean(x_delta_vec);
[~,DT_Ind] = sortrows(DataTable(:, 1:4));
DataTable = DataTable(DT_Ind,:);

if ~isempty(DS)
    DataTypes = DataTable.('Data Type');
    BoolVec = false(size(DataTypes));
    if iscell(DS)
        for c = DS
            BoolVec = BoolVec | DataTypes==c{:};
        end
    elseif ChrStr(DS)
        BoolVec = DataTypes==DS;
    end
    DataTable = DataTable(BoolVec,:);
end
if TeRflag
    DataTable = DataTable(DataTable.T>=TR(1) & DataTable.T<=TR(2),:);
else
    DataTable = DataTable(DataTable.T<TR(1) | DataTable.T>TR(2),:);
end
if ~isempty(SN)
    DataTable = DataTable(ismember(DataTable.('Set Number'),SN),:);
end
end

