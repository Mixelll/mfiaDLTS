function [CtCellsFrames, tDelta, TimeTemperature] = MFIA_DAQ_MatchData(Ct_path, TimeZone, StartTime, EndTime, DataSelect)
DataSelectC = fn_struct2cell(DataSelect);
DataSelectC(5,:) = {'Voltage', 'Resistance', 'Capacitance'};
CtCellsFrames = cell(3,0);
xls = dir([Ct_path '\*.xls']);
if length(xls)==1
    LogT = readmatrix([Ct_path '\' xls.name]);
    LogT_time = ParseDate(readcell([Ct_path '\' xls.name], 'Range', [2 2 2 2]))*1e6;
elseif isempty(xls)
    error('Temperature log file (.xls) not found')
else 
    error('More than 1 Temperature log files (.xls) were found')
end
mat = dir([Ct_path '\*.mat']);
if length(xls)==1
    File = open([Ct_path '\' mat.name]);
elseif isempty(xls)
    error('MFIA data log file (.mat) not found')
else 
    error('More than 1 MFIA data log files (.mat) were found')
end
TimeTemperature = [1000*LogT(:,1)+LogT_time LogT(:,5)];
for c = DataSelectC
    if c{2}
        SelectedData.(c{5}) = eval(['File.dev5168' c{3}]);
    end
end
SD_Fnames = fieldnames(SelectedData)';
N_Frames = length(SelectedData.(SD_Fnames{1}));
for i=1:N_Frames
    CurveStruct = SelectedData.(SD_Fnames{1}){i};
    DataC = CurveStruct.value;
    tDelta = CurveStruct.header.gridcoldelta;
    tOffset = CurveStruct.header.gridcoloffset;
    if tOffset<0
        DataC = DataC(ceil(-tOffset/tDelta)+1:end);
    end
    TransposeFlag = 0;
    if size(DataC,2)>size(DataC,1)
        TransposeFlag = 1;
        DataC = DataC.';
    end
    DataTime = StartTime+tDelta*(1:length(DataC))';
    DataC_time = [DataTime(DataTime<EndTime) DataC(DataTime<EndTime)];
    Time = double(CurveStruct.header.systemtime) + TimeZone*3600*1e6;
    T_Curve = TimeTemperature(find((TimeTemperature(:,1)-Time)>0,1),2);
    DataStruct.(SD_Fnames{1}) = DataC_time;
    for c = SD_Fnames(2:end)
        DataC = SelectedData.(c{:}){i}.value(ceil(-tOffset/tDelta)+1:end);
        if TransposeFlag
            DataC = DataC.';
        end
        DataStruct.(c{:}) = [DataTime(DataTime<EndTime) DataC(DataTime<EndTime)];
    end
    CtCellsFrames(:,end+1) = {['T=' num2str(T_Curve,4)] ; DataStruct ; {'T' T_Curve}};
end
TimeTemperature(:,1) = (TimeTemperature(:,1) - LogT_time)/1e6;
end

