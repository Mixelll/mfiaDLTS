function OutStruct = MFIA_maxmin_axes(InStruct)
OutStruct.frequency = [nan nan];
OutStruct.amplitude = [nan nan];
OutStruct.offset = [nan nan];
OutStruct.log_frequency = [nan nan];
if any(strcmpi(fieldnames(InStruct), 'axes'))
    for s = InStruct
        OutStruct.frequency = [min(OutStruct.frequency(1),min(s.axes.frequency)) max(OutStruct.frequency(2),max(s.axes.frequency))];
        OutStruct.amplitude = [min(OutStruct.amplitude(1),min(s.axes.amplitude)) max(OutStruct.amplitude(2),max(s.axes.amplitude))];
        OutStruct.offset = [min(OutStruct.offset(1),min(s.axes.offset)) max(OutStruct.offset(2),max(s.axes.offset))];
        OutStruct.frequency = [min(OutStruct.log_frequency(1),min(s.axes.log_frequency)) max(OutStruct.log_frequency(2),max(s.axes.log_frequency))];
    end
else
        for s = InStruct
        OutStruct.frequency = [min(OutStruct.frequency(1),min(s.frequency)) max(OutStruct.frequency(2),max(s.frequency))];
        OutStruct.amplitude = [min(OutStruct.amplitude(1),min(s.amplitude)) max(OutStruct.amplitude(2),max(s.amplitude))];
        OutStruct.offset = [min(OutStruct.offset(1),min(s.offset)) max(OutStruct.offset(2),max(s.offset))];
        OutStruct.frequency = [min(OutStruct.log_frequency(1),min(s.log_frequency)) max(OutStruct.log_frequency(2),max(s.log_frequency))];
        end
end

