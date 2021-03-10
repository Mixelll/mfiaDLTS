function [DataTableSegments, DataStructCells, x_delta] = MFIA_DAQ_Process(DataTable, varargin) 
CellStrChar = @(s) iscell(s) || isstring(s) || ischar(s);
ChrStr = @(s) ischar(s) || isstring(s);
CellStrCharFn = @(s) iscell(s) || isstring(s) || ischar(s) || isa(s, 'function_handle');
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('DataSelect', {'Capacitance'}, CellStrChar); % Select data type
p.addParameter('TemperatureRange', [-inf inf 1], @isnumeric); % Select temperature range. Third element: true to select inside range, false for outside 
p.addParameter('SetNumbers', [], @isnumeric); % Select sets by numbers
p.addParameter('Range', [-inf inf], @isnumeric); % Trim data by time range
p.addParameter('CurveSelect', {}, CellStrCharFn); % Select sets by name['F' PulseVoltage 'T' MeasureVoltage]
p.addParameter('MovMean', 1, @isnumeric); % Apply a moving average
p.addParameter('DeleteFalseTransition', true); % Replace data with NaN after a false transition

p.parse(varargin{:});

DS = p.Results.DataSelect;
SN = p.Results.SetNumbers;
R = p.Results.Range;
CS = p.Results.CurveSelect;
MM = p.Results.MovMean;
TR = p.Results.TemperatureRange;
DFT = p.Results.DeleteFalseTransition;
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

DSCflag = true;
if length(DS)~=1
    warning(['DataSelect has ' num2str(length(DS)) ' entries, DataStructCells will not be returned'])
    DSCflag = false;
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
if TeRflag
    DataTable = DataTable(DataTable.T>=TR(1) & DataTable.T<=TR(2),:);
else
    DataTable = DataTable(DataTable.T<TR(1) | DataTable.T>TR(2),:);
end
if ~isempty(SN)
    DataTable = DataTable(ismember(DataTable.('Set Number'),SN),:);
end

SetNames = DataTable.('Set Name').';
SetNumbers = unique(DataTable.('Set Number')).';
Curves = {};
Funcs = {};
x_deltaArray = [];
DataStructCells = struct;
DataTableSegments = table('Size',[0,14],'VariableNames', {'Data Type', 'Set Number', 'Set Name', 'T', 'FromTo','FromToVar', 't0', 't end', 'Valid Data Points', 'MovMean', 'Mean', 'Var', 'Norm SD', 'Data'},...
    'VariableTypes', {'string','double','string','double','string','string','double','double','double','double','double','double','double','cell'});
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
                FalseTrans = [];
                for s = Curves{i}
                    try 
                        FnInput(end+1) = DT(DT.('FromToVar')==s,:).Data;
                        FalseTrans = [FalseTrans DT(DT.('FromToVar')==s,:).('False Transitions'){:}];
                    catch 
                        error(['Expected one output from table but instead got ' num2str(height(DT(DT.('FromToVar')==s,:).Data)) ' outputs'])
                    end
                end
            end
            MinL = min(cellfun(@(c)length(c), FnInput));
            x = FnInput{1}(1:MinL,1); x_delta =  mean(diff(x)); x_deltaArray(end+1) = x_delta;
            FnInput = cellfun(@(c) c(1:MinL,2), FnInput , 'UniformOutput',false);
            FnOutput = Funcs{i}(FnInput{:});
            if DFT(1)
                if ~isempty(FalseTrans)
                    MFT = min(FalseTrans);
                    FnOutput(MFT:end) = [];
                    x(MFT:end) = [];
                end
            end
            Rvec = ~excludedata(x, FnOutput, 'box',R); x = x(Rvec); FnOutput = FnOutput(Rvec);
            NanVec = isnan(FnOutput);
            VDP = length(x) - sum(NanVec); % Valid data points
            FnOutput = movmean(FnOutput,MM,'omitnan');
            Data = {[x FnOutput]};
            FnStr = func2str(Funcs{i});
            FromToVar = regexprep(FnStr(find(FnStr==')',1)+1:end), {'-'}, {'M'});
            FromTo = regexprep(FromToVar, {'m','M','p'}, {'-','-','\.'});
            SetName = unique(DT.('Set Name'));
            T = unique(DT.T);
            Mean = mean(FnOutput(~NanVec));
            Var = moment(FnOutput(~NanVec),2);
            NSD = sqrt(Var)/Mean;
            DataTableSegments(end+1,:) = {d set SetName T FromTo FromToVar x(1) x(end) VDP MM Mean Var NSD Data};
            if DSCflag
                SetNameChar = convertStringsToChars(SetName);
                if isfield(DataStructCells, FromToVar)
                    DataStructCells.(FromToVar)(:,end+1) = {SetNameChar; Data{:}; {'T' T}};
                else
                    DataStructCells.(FromToVar) = {SetNameChar; Data{:}; {'T' T}};
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

