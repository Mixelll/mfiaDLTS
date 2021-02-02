function [DataStruct] = MFIA_DAQ_FrameSegmenter(DataCells,varargin)
p = inputParser;
p.KeepUnmatched=true;
p.addParameter('TransitionTime', 1.8667e-04, @isnumeric); % Actual voltage transition time should be less than TransitionTime
p.addParameter('RoundVoltage', 5e-2, @isnumeric); % Round displayed voltage to within RoundVoltage
p.addParameter('VoltageTransitionCriteria', 1e-2, @isnumeric); % Changes of voltage larger than VoltageTransitionCriteria will be treated as transition points
p.addParameter('UniformInput', 0.1, @isnumeric); % Discard transition points if they appear less frequently than UniformInput fraction total points
p.addParameter('CoerceVoltage', 0.1, @isnumeric); % If before-segment-voltage (Vfrom) or segment-voltage (Vto) of a certain transition point is within 
% CoerceVoltage from the expected voltage (mean - from the majority of data sets), it will be coerced to the expected voltage
p.parse(varargin{:});

Ncells = size(DataCells,2);
Fnames = fieldnames(DataCells{2,1})';
VfieldName = Fnames{contains(Fnames, 'volt', 'IgnoreCase',true)};
RV = p.Results.RoundVoltage; % Voltage round bracket
TC = p.Results.VoltageTransitionCriteria; % Transition Criteria
UI = p.Results.UniformInput;
CV = p.Results.CoerceVoltage;

SetTable = table('Size',[0,7],'VariableNames', {'Set Number', 'Set Name', 'Transition Index', 'From', 'To', 'From Round', 'To Round'},...
    'VariableTypes', {'double','string','double','double','double','double','double'});
TransitionTable = table('Size',[0,7],'VariableNames', {'Transition Index', 'Occurances', 'Set Names', 'From Round', 'From SD', 'To Round', 'To SD'},...
    'VariableTypes', {'double','double','cell','double','double','double','double'});
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
    TransitionTable(end+1,:) = {i Nsets {SetNames} CustomRound(mean(FromArray),RV) std(FromArray) CustomRound(mean(ToArray),RV) std(ToArray)};
end
disp(TransitionTable);

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
        if ~UI || sum(TransitionTable(TransitionTable.('Transition Index')==i,:).Occurances/Ncells)>=UI
            if isempty(Tvec) || i<(Tvec(end)+TP)
                Tvec(end+1) = i;
            end
            if ~isempty(Tvec) && (i>(Tvec(end)+TP) || i==yVdiffInd(end))
                TransitionIndex = Tvec(ceil(length(Tvec)/2));
                Vfrom = mean(yV(max(Tvec(1)-2*TP+1,1):Tvec(1)-TP)); VfromExpected = TransitionTable(TransitionTable.('Transition Index')==TransitionIndex,:).('From Round');
                if abs(Vfrom-VfromExpected)<CV, Vfrom = VfromExpected; end
                Vto = mean(yV(Tvec(end)+TP:min(Tvec(end)+2*TP-1,SL))); VtoExpected = TransitionTable(TransitionTable.('Transition Index')==TransitionIndex,:).('From Round');
                if abs(Vto-VtoExpected)<CV, Vto = VtoExpected; end
                
                if ~isempty(Tcell)
                    Tcell{3,end} = {[Tcell{3,end} TransitionIndex]};
                end
                Tcell(:,end+1) = {['from ' num2str(CustomRound(Vfrom,RV)) ' to ' num2str(CustomRound(Vto,RV))] ; [Vfrom Vto] ; TransitionIndex+1};
                Tvec = i;    
            end
        end
    end
    if ~isempty(Tcell)
        Vfrom1 = Tcell{2,end}(2);
        Vto1 = Tcell{2,1}(1);
        Tcell = [{['from ' num2str(CustomRound(Vfrom1,RV)) ' to ' num2str(CustomRound(Vto1,RV))] ; [Vfrom1 Vto1]; [1 Tcell{3,1}(1)-1]} Tcell];
        Tcell(3,end) = {[Tcell{3,end} SL]};
    else
        Vto = mean(yV(end-TP:end-1));
        if ~isempty(varargin)
            if length(varargin{:})==2
                [~, FromInd] = max(abs(varargin{:}-Vto));
                Vfrom = varargin{:}(FromInd);
            elseif length(varargin{:})==1
                Vfrom = varargin{:};
            end
        else
            Vfrom = 0;
        end
        Tcell = {['from ' num2str(CustomRound(Vfrom,RV)) 'to ' num2str(CustomRound(Vto,RV))] ; [Vfrom Vto] ; [1 SL]};
    end
    
    DataStruct = struct;
    for c = Fnames
        DataC = Data.(c{:});
        for tc = Tcell
            FT = tc{2}; % From V to V
            ST = tc{3}; % Segment start-stop
            DataC_Seg = DataC(ST(1):ST(2),2); xDataC_Seg = x_delta*(1:length(DataC_Seg)).';
            TempStruct.(regexprep(['F' num2str(CustomRound(FT(1),RV)) 'T' num2str(CustomRound(FT(2),RV),3)], {'-','\.'}, {'m','p'})) = [xDataC_Seg DataC_Seg];
        end
        DataStruct.(c{:}) = TempStruct;
    end
 
end
end

