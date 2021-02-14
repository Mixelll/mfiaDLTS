function [DataTable, DataCellsStruct, TransTable, SetTable] = MFIA_DAQ_FrameSegmenter(DataCells, varargin)
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('TransitionTime', 1.8667e-04, @isnumeric); % Actual voltage transition time should be less than TransitionTime
p.addParameter('RoundVoltage', 5e-2, @isnumeric); % Round displayed voltage to within RoundVoltage
p.addParameter('VoltageTransitionCriteria', 1e-2, @isnumeric); % Changes of voltage larger than VoltageTransitionCriteria will be treated as transition points
p.addParameter('UniformInput', 0.1, @isnumeric); % Discard transition points if they appear less frequently than UniformInput fraction total points
p.addParameter('VoltageValues', [], @isnumeric); % Input expected voltage values. The function will try to assign a missing value in case of partially captured frames
p.addParameter('CoerceVoltage', 0.12, @isnumeric); % If before-segment-voltage (Vfrom) or segment-voltage (Vto) of a certain transition point is within 
% CoerceVoltage from the expected voltage (mean - from the majority of data sets), it will be coerced to the expected voltage
p.parse(varargin{:});

Ncells = size(DataCells,2);
Fnames = fieldnames(DataCells{2,1})';
VfieldName = Fnames{contains(Fnames, 'volt', 'IgnoreCase',true)};
RV = p.Results.RoundVoltage; % Voltage round bracket
TC = p.Results.VoltageTransitionCriteria; % Transition Criteria
UI = p.Results.UniformInput;
CV = p.Results.CoerceVoltage;
VV = p.Results.VoltageValues;

DataCellsStruct = DataCells;
DataTable = table('Size',[0,12],'VariableNames', {'Data Type', 'Set Number', 'Set Name', 'T', 'FromTo','FromToVar', 'From', 'To', 't0', 'Length Seconds', 'Length', 'Mean', 'Var', 'Norm SD', 'Data'},...
    'VariableTypes', {'string','double','string','double','string','string','double','double','double','double','double','double','double','double','cell'});
SetTable = table('Size',[0,7],'VariableNames', {'Set Number', 'Set Name', 'Transition Index', 'From', 'To', 'From Round', 'To Round'},...
    'VariableTypes', {'double','string','double','double','double','double','double'});
TransTable = table('Size',[0,7],'VariableNames', {'Transition Index', 'Occurances', 'Set Names', 'From Round', 'From SD', 'To Round', 'To SD'},...
    'VariableTypes', {'double','double','cellstr','double','double','double','double'});
for ic = 1:Ncells
    Data = DataCells{2,ic};
    DataV = Data.(VfieldName);
    xV = DataV(:,1); yV = DataV(:,2); x_delta = mean(diff(xV));

    TP = round(p.Results.TransitionTime/x_delta); % Transition points length
    SL = length(yV);
    yVdiff = diff(yV); yVdiffInd = find(abs(yVdiff)>TC).'+1; yVdiffInd(yVdiffInd<TP) = []; yVdiffInd(yVdiffInd>(SL-TP)) = []; % Find transition indexes

    for i = yVdiffInd
        Vfrom = mean(yV(max(i-2*TP+1,1):i-TP)); Vto = mean(yV(i+TP:min(i+2*TP-1,SL)));
        SetTable(end+1,:) = {ic DataCells{1,ic} i Vfrom Vto CustomRound(Vfrom,RV) CustomRound(Vto,RV)};
    end
end
for i = unique(SetTable.('Transition Index')).'
    TempTable = SetTable(SetTable.('Transition Index')==i,:);
    Nsets = height(TempTable);
    SetNames = TempTable.('Set Name');
    FromArray = TempTable.('From');
    ToArray = TempTable.('To');
    TransTable(end+1,:) = {i Nsets {SetNames} CustomRound(mean(FromArray),RV) std(FromArray,1) CustomRound(mean(ToArray),RV) std(ToArray,1)};
end
[~,DT_Ind] = sortrows(DataTable(:,1:7));
DataTable = DataTable(DT_Ind,:);
disp(TransTable);
DisInd = [];
DisSet = {};
NoTransVtoVec = [];
for ic = 1:Ncells
    Data = DataCells{2,ic};
    DataV = Data.(VfieldName);
    xV = DataV(:,1); yV = DataV(:,2); x_delta = mean(diff(xV));
    
    TP = round(p.Results.TransitionTime/x_delta); % Transition points length
    SL = length(yV);
    yVdiff = diff(yV); yVdiffInd = find(abs(yVdiff)>TC).'+1; yVdiffInd(yVdiffInd<TP) = []; yVdiffInd(yVdiffInd>(SL-TP)) = []; % Find transition indexes
    
    if isempty(yVdiffInd), yVdiffInd = []; end
    Tvec = []; Tcell = cell(3,0); % Initialize temporary containers
    for i = yVdiffInd
        if ~UI || sum(TransTable(TransTable.('Transition Index')==i,:).Occurances/Ncells)>=UI
            if isempty(Tvec) || i<(Tvec(end)+TP)
                Tvec(end+1) = i;
            end
            if ~isempty(Tvec) && (i>(Tvec(end)+TP) || i==yVdiffInd(end))
                TransInd= Tvec(ceil(length(Tvec)/2));
                Vfrom = mean(yV(max(Tvec(1)-2*TP+1,1):Tvec(1)-TP)); VfromC = Vfrom; VfromExp = TransTable(TransTable.('Transition Index')==TransInd,:).('From Round');
                if abs(Vfrom-VfromExp)<=CV, VfromC = VfromExp; end
                Vto = mean(yV(Tvec(end)+TP:min(Tvec(end)+2*TP-1,SL))); VtoC = Vto; VtoExp = TransTable(TransTable.('Transition Index')==TransInd,:).('To Round');
                if abs(Vto-VtoExp)<=CV, VtoC = VtoExp; end
                
                if ~isempty(Tcell)
                    Tcell{3,end} = {[Tcell{3,end} TransInd+1]};
                end
                Tcell(:,end+1) = {[VfromC VtoC] ; [Vfrom Vto] ; TransInd+1};
                Tvec = i;    
            end
        else
            DisInd(end+1) = i;
            DisSet(end+1) = DataCells(1,ic);
            
        end
    end
    if ~isempty(Tcell)
        Vfrom1 = Tcell{2,end}(2);
        Vto1 = Tcell{2,1}(1);
        Vfrom1C = Tcell{1,end}(2);
        Vto1C = Tcell{1,1}(1);
        Tcell = [{[Vfrom1C Vto1C]; [Vfrom1 Vto1]; [1 Tcell{3,1}(1)-1]} Tcell];
        Tcell(3,end) = {[Tcell{3,end} SL]};
    else
        Vto = mean(yV(end-TP:end-1)); VtoC = Vto; NoTransVtoVec(end+1) = Vto; NoTransVtoVecAvg = mean(NoTransVtoVec);
        if abs(Vto-NoTransVtoVecAvg)<CV, VtoC = CustomRound(NoTransVtoVecAvg,RV); end
        if ~isempty(VV)
            if length(VV)==2
                [~, FromInd] = max(abs(VV-Vto));
                Vfrom = VV(FromInd); VfromC = Vfrom;  VVt = VV; VVt(FromInd) = []; VtoC = VVt;
            elseif length(VV)==1
                Vfrom = VV;
            end
        else
            Vfrom = 0;
        end
        Tcell = {[VfromC VtoC] ; [Vfrom Vto] ; [1 SL]};
    end
    
    DataStruct = struct;
    for c = Fnames
        DataC = Data.(c{:});
        for tc = Tcell
            FTC = tc{1}; % From V to V coerced
            FT = tc{2}; % From V to V actual
            ST = tc{3}; % Segment start-stop
            DataC_Seg = DataC(ST(1):ST(2),2); xDataC_Seg = x_delta*(1:length(DataC_Seg)).';
            FTstr = ['F' num2str(CustomRound(FTC(1),RV)) 'T' num2str(CustomRound(FTC(2),RV),3)]; FTvar = regexprep(FTstr, {'-','\.'}, {'m','p'});
            TempStruct.(FTvar) = [xDataC_Seg DataC_Seg];
            Mean = mean(DataC_Seg);
            Var = moment(DataC_Seg,2);
            NSD = sqrt(Var)/Mean;
            DataTable(end+1,:) = {c{:} ic DataCells{1,ic} DataCells{3,ic}{2} FTstr FTvar CustomRound(FT(1),RV) CustomRound(FT(2),RV) xDataC_Seg(1) xDataC_Seg(end) length(xDataC_Seg) Mean Var NSD {[xDataC_Seg DataC_Seg]}};
        end
        DataStruct.(c{:}) = TempStruct;
    end
 DataCellsStruct{2,ic} = DataStruct; 
end
[~,DT_Ind] = sortrows(DataTable(:,1:7));
DataTable = DataTable(DT_Ind,:);
if ~isempty(DisInd)
    DisSet = unique(DisSet);
    DisSetStr = [];
    for c = DisSet
        DisSetStr = [DisSetStr ', ' c{:}];
    end 
    warning(['Discarded a total of ' num2str(length(DisInd)) ' false transitions across ' num2str(length(DisSet)) ' sets' ...
        newline 'Discarded transition times: ' sprintf('%0.1g ',unique(DisInd)*x_delta) '[sec]'...
        newline 'From datasets: ' DisSetStr]);
end
end

