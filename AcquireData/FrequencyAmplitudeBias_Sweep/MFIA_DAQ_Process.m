function [DataTableOut, DataStructCells, x_delta] = MFIA_DAQ_Process(DataTable, varargin) 
CellStrChar = @(s) iscell(s) || isstring(s) || ischar(s);
ChrStr = @(s) ischar(s) || isstring(s);
CellStrCharFn = @(s) iscell(s) || isstring(s) || ischar(s) || isa(s, 'function_handle');
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('DataSelect', {'Capacitance'}, CellStrChar); %
p.addParameter('TemperatureRange', [-inf inf], @isnumeric); %
p.addParameter('NumberRange', [1 inf], @isnumeric); %
p.addParameter('TimeRange', [-inf inf], @isnumeric); % 
p.addParameter('CurveSelect', {}, CellStrCharFn); % 
p.addParameter('MovMean', 1, @isnumeric); % 
 
p.parse(varargin{:});

DS = p.Results.DataSelect;
TR = p.Results.TimeRange;
NR = p.Results.NumberRange;
CS = p.Results.CurveSelect;
MM = p.Results.MovMean;
DSCflag = 1;
if length(DS)~=1
    warning(['DataSelect has ' num2str(length(DS)) ' entries, DataStructCells will not be returned'])
end
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
DataTable = DataTable(DataTable.T>=TR(1) & DataTable.T<=TR(2),:);
DataTable = DataTable(DataTable.('Set Number')>=NR(1) & DataTable.('Set Number')<=NR(2),:);
SetNames = DataTable.('Set Name').';
SetNumbers = unique(DataTable.('Set Number')).';
Curves = {};
Funcs = {};
x_deltaArray = [];
DataStructCells = struct;
DataTableOut = table('Size',[0,10],'VariableNames', {'Data Type', 'Set Number', 'Set Name', 'T', 'FromTo','FromToVar', 't0', 'Length Seconds','MovMean', 'Data'},...
    'VariableTypes', {'string','double','string','double','string','string','double','double','double','cell'});
if ~isempty(CS)
    if iscell(CS)
        for c = CS
            cst = c{:};
            if ChrStr(cst)
            cst = regexprep(cst, {'-','\.'}, {'m','p'});
            Curves{end+1} = convertCharsToStrings(cst);
            CSchar = convertStringsToChars(cst);
            Funcs{end+1} = str2func(['@(' CSchar ')' CSchar]);
            elseif isa(cst, 'function_handle')
                Funcs{end+1} = cst;
                FnStr = func2str(cst);
                I1 = find(FnStr=='(',1);
                I2 = find(FnStr==')',1);
                Curves{end+1} = cellfun(@(c)convertCharsToStrings(c),strsplit(FnStr(I1+1:I2-1),','));
            end
        end
    elseif ChrStr(CS)
            CS = regexprep(CS, {'-','\.'}, {'m','p'});
            Curves{end+1} = convertCharsToStrings(CS);
            CSchar = convertStringsToChars(CS);
            Funcs{end+1} = str2func(['@(' CSchar ')' CSchar]);
    elseif isa(CS, 'function_handle')
            Funcs{end+1} = CS;
            FnStr = func2str(CS);
            I1 = find(FnStr=='(',1);
            I2 = find(FnStr==')',1);
            Curves{end+1} = cellfun(@(c)convertCharsToStrings(c),strsplit(FnStr(I1+1:I2-1),','));
    end
    for d = unique(DataTable.('Data Type')).' 
        DataTableD = DataTable(DataTable.('Data Type')==d,:);
        for set = SetNumbers
            DT = DataTableD(DataTableD.('Set Number')==set,:);
            for i = 1:length(Curves)
                FnInput = {};
                for s = Curves{i}
                    try 
                        FnInput(end+1) = DT(DT.('FromToVar')==s,:).Data;
                    catch 
                        error(['Expected one output from table but instead got ' num2str(height(DT(DT.('FromToVar')==s,:).Data)) ' outputs'])
                    end
                end
            end
            MinL = min(cellfun(@(c)length(c), FnInput));
            x = FnInput{1}(1:MinL,1); x_delta =  mean(diff(x)); x_deltaArray(end+1) = x_delta;
            FnInput = cellfun(@(c) c(1:MinL,2), FnInput , 'UniformOutput',false);
            FnOutput = movmean(Funcs{i}(FnInput{:}),MM);
            Data = {[x FnOutput]};
            FnStr = func2str(Funcs{i});
            FromToVar = FnStr(find(FnStr==')',1)+1:end);
            FromTo = regexprep(FromToVar, {'m','p'}, {'-','\.'});
            SetName = unique(DT.('Set Name'));
            T = unique(DT.T);
            DataTableOut(end+1,:) = {d set SetName T FromTo FromToVar x(1) x(end) MM Data};
            if DSCflag
                if isfield(DataStructCells, FromToVar)
                    DataStructCells.(FromToVar)(:,end+1) = {SetName; Data; {'T' T}};
                else
                    DataStructCells.(FromToVar) = {SetName; Data; {'T' T}};
                end
            end
        end
    end
x_delta = mean(x_deltaArray);   
end
        
%             if ischar(cs)
%                 cs = convertCharsToStrings(cs);
%             end
%             for s = cs
%                 if strcmpi(s,{'-', '+', '*', '/', '.^'})
%                 else
                    
                     
end

