function [Cells, CellsofStruct,x_delta] = MFIA_average_t(DataCells,N_signals,SignalSlect,SegmentSelect,MovMean)
if isempty(MovMean)
    MovMean = 1;
end
CellsofStruct = DataCells;
Cells = {};
for ic = 1:size(DataCells,2)
    Data = DataCells{2,ic};
    SignalLength = length(Data)/N_signals;
    SignalCells = mat2cell(Data,repmat(SignalLength,1,N_signals),2);
    Signals = {};
    for sc = SignalCells'
        DataOne = sc{:};
        x = DataOne(:,1);
        y = DataOne(:,2);
        if sum(abs(y)>=11)>SignalLength/2
            Signals{end+1} = 'Resistance';
        elseif sum(abs(y)<=0.000001)>SignalLength/2
            Signals{end+1} = 'Capacitance';
        elseif sum(abs(y)>=0.00001 & abs(y)<=10)>SignalLength/4
            xV = x;
            yV = y;
            Signals{end+1} = 'Voltage';
        else
            Signals{end+1} = 'Unspecified';
        end
    end
    x_delta = mean(diff(xV));
    yMean = mean(y);
    yLogicalHigh = yV>yMean;
%     yLogicalLow = yV<yMean;
    Border = [0 ; abs(diff(yLogicalHigh))];
    BorderFind = find(Border);
    Start = BorderFind(1);
    End = BorderFind(end)-1;
    N_segments = sum(Border)-1;
    SL = floor(mean(diff(BorderFind)));
    
    strFirstV = num2str(mean(yV(BorderFind(1)+1:BorderFind(2)-1)));
    strSecondV = num2str(mean(yV(BorderFind(2)+1:BorderFind(3)-1)));
    
    HL_Criteria = all(yLogicalHigh(BorderFind(1)+1:BorderFind(2)-1));
    if HL_Criteria
        StrFirst = 'L2High';
        StrSecond = 'H2Low';
    else
        StrFirst = 'H2Low';
        StrSecond = 'L2High';
    end
    D = struct;
    for isc = 1:length(SignalCells)
        if any(strcmpi(Signals{isc},SignalSlect))
            DataOne = SignalCells{isc};
            y = DataOne(:,2);
            yFirst = zeros(SL,ceil(N_segments/2));
            ySecond = zeros(SL,floor(N_segments/2));
            i = 1;
            for i_Border = BorderFind(1:2:end)'
                if i_Border+SL>length(y)
                    break
                end
                yFirst(:,i) = y(i_Border:i_Border+SL-1);
                i = i+1;
            end
            i = 1;
            for i_Border = BorderFind(2:2:end)'
                if i_Border+SL>length(y)
                    break
                end
                ySecond(:,i) = y(i_Border:i_Border+SL-1);
                i = i+1;
            end
            yFirstMean = movmean(mean(yFirst,2),MovMean);
            ySecondMean = movmean(mean(ySecond,2),MovMean);
            xVec = (0:x_delta:x_delta*(SL-1)).';
            
            First = [xVec yFirstMean];
            Second = [xVec ySecondMean];
            if contains(StrFirst,SegmentSelect,'IgnoreCase',true)
                D.(Signals{isc}).(StrFirst) = First;
                Cells(:,end+1) = [{[DataCells{1,ic} ' ' Signals{isc} ' ' StrFirst ' V=' strFirstV]}; {First}; DataCells(3,ic)];
            end
            if contains(StrSecond,SegmentSelect,'IgnoreCase',true)
                D.(Signals{isc}).(StrSecond) = Second;
                Cells(:,end+1) = [{[DataCells{1,ic} ' ' Signals{isc} ' ' StrSecond ' V=' strSecondV]}; {Second}; DataCells(3,ic)];
            end
        end
    end
    CellsofStruct{2,ic} = D;  
end
end
