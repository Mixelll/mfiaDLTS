function [select_data, full_data, OutSettings]  = MFIA_freq_amp_bias_sweep(device_id, additional_settings, SetByRange, sweep_order, sweep_range, pts, frequency_vec, amplitude_vec, offset_vec, read_param_struct, varargin)
% Perform a frequency/test_signal (amplitude)/bias_voltage (offset) sweep and gather demodulator data.
%
% NOTE Please ensure that the ziDAQ folders 'Driver' and 'Utils' are in your
% Matlab path. To do this (temporarily) for one Matlab session please navigate
% to the ziDAQ base folder containing the 'Driver', 'Examples' and 'Utils'
% subfolders and run the Matlab function ziAddPath().
% >>> ziAddPath;
%
% Use either of the commands:
% >>> help ziDAQ
% >>> doc ziDAQ
% in the Matlab command window to obtain help on all available ziDAQ commands.
vararginTOP = varargin;
lf = length(frequency_vec);
la = length(amplitude_vec);
lo = length(offset_vec);
lpairs = eval(['max(length(' sweep_order{2} '_vec),length(' sweep_order{3} '_vec))']);
if ~(sum([lf la lo]>0)>=2 && (lf==la || lf==lo || lo==la))
    fprintf('sweep_order #2 and #3 vecs do not match in length');
    return
end
SetByRange(5,:) = {true};   
clear ziDAQ;

if ~exist('device_id', 'var')
    error(['No value for device_id specified. The first argument to the ' ...
           'example should be the device ID on which to run the example, ' ...
           'e.g. ''dev2006'' or ''uhf-dev2006''.'])
end

% Check the ziDAQ MEX (DLL) and Utility functions can be found in Matlab's path.
if ~(exist('ziDAQ', 'file') == 3) && ~(exist('ziCreateAPISession', 'file') == 2)
    fprintf('Failed to either find the ziDAQ mex file or ziDevices() utility.\n')
    fprintf('Please configure your path using the ziDAQ function ziAddPath().\n')
    fprintf('This can be found in the API subfolder of your LabOne installation.\n');
    fprintf('On Windows this is typically:\n');
    fprintf('C:\\Program Files\\Zurich Instruments\\LabOne\\API\\MATLAB2012\\\n');
    return
end

% The API level supported by this example.
apilevel_example = 6;
% Create an API session; connect to the correct Data Server for the device.
[device, props] = ziCreateAPISession(device_id, apilevel_example);
ziApiServerVersionCheck();

branches = ziDAQ('listNodes', ['/' device ], 0);
if ~any(strcmpi(branches, 'DEMODS'))
  fprintf('\nThis example requires lock-in functionality which is not available on %s.\n', device);
  return
end

% Define parameters relevant to this example. Default values specified by the
% inputParser below are overwritten if specified as name-value pairs via the
% `varargin` input argument.
p = inputParser;
p.KeepUnmatched=true;
isnumscalar = @(x) isnumeric(x) && isscalar(x);
isnonnegscalar = @(x) isnumeric(x) && isscalar(x) && (x > 0);

% IA precision -> measurement speed: 0 - low->fast, 1 - high->medium,
% 2 - very high->slow
p.addParameter('IA_precision', 1, isnumscalar);

% Set IA parameter extraction model (from impedance). 0 - Rp Cp,
% 1 - Rs Cs, 2 - Rs Ls, 3 - G B, 4 -D Cs, 5 - Qs Cs, 6 - D Ls,
% 7 - Q Ls, 8 - Rp Lp, 9 - D Cp
p.addParameter('model', 0, isnumscalar);

% Enable two terminal measurement.
p.addParameter('two_terminal', 1, isnumscalar);

% Enable two terminal cable length compensation [m].
p.addParameter('cable_length', 1, isnumscalar);

% Enable high-pass filter.
p.addParameter('AC', 0, isnumscalar);

% Enable four terminal high-pass filter.
p.addParameter('imp50ohm', 0, isnumscalar);

% Enable current auto range.
p.addParameter('auto_range', 0, isnumscalar);

% Current manual range, [A].
p.addParameter('current_range', 10e-6, @isnumeric);

% Voltage manual range, [V].
p.addParameter('voltage_range', 3, @isnumeric);

% demod time constant, [s].
p.addParameter('demod_time_constant', 0.007, @isnumeric);

% demod data transfer rate, [Sa/sec].
p.addParameter('demod_rate', 13.39e3, @isnumeric);

% For IA:  Selects the filter roll off to use for the sweep in fixed bandwidth mode. Range between 6 dB/oct and 48 dB/ oct.
p.addParameter('demod_LFP_order', 8, isnonnegscalar);

p.parse(varargin{:});
unmatched_vars = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)];
unmatched_vars = unmatched_vars.';
OutSettings.num = update_structure(p.Results,p.Unmatched);
%% Set default additional settings
additional_settings_internal.enable_default = true;
% Enable DLCP condition on offset (Vbias): ActualOffset = InputOffset - Amplitude/2
additional_settings_internal.output.DLCP = false;
% Define device channels.
additional_settings_internal.channels.demod_c = '0'; % demod channel, for paths on the device
additional_settings_internal.channels.demod_idx = str2double(additional_settings_internal.channels.demod_c)+1; % 1-based indexing, to access the data
additional_settings_internal.channels.out_c = '0'; % signal output channel
% Get the value of the instrument's default Signal Output mixer channel.
additional_settings_internal.channels.out_mixer_c = num2str(ziGetDefaultSigoutMixerChannel(props, str2double(additional_settings_internal.channels.out_c)));
additional_settings_internal.channels.in_c = '0'; % signal input channel
additional_settings_internal.channels.osc_c = '0'; % oscillator
additional_settings_internal.channels.imp_c = '0'; % IA channel
additional_settings_internal.channels.imp_index = str2double(additional_settings_internal.channels.imp_c)+1; % IA, 1-based indexing, to access the data
% Graphs
additional_settings_internal.display.graph.disp = true;
additional_settings_internal.display.graph.during_sweep = false;
additional_settings_internal.display.graph.save.if = false;
additional_settings_internal.display.graph.save.path = '';
% Text
additional_settings_internal.display.text.major.disp = true;
additional_settings_internal.display.text.major.each_sweep = true;
additional_settings_internal.display.text.minor.disp = true;
additional_settings_internal.display.text.minor.each_sweep = false;

% Overwrite default additional settings
additional_settings_internal = update_structure(additional_settings_internal, additional_settings); 

c = additional_settings_internal.channels;

major.disp = additional_settings_internal.display.text.major.disp;
minor.disp = additional_settings_internal.display.text.minor.disp;
enable_default = additional_settings_internal.enable_default;
OutSettings.additional = additional_settings_internal;
% Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
ziDisableEverything(device);

%% Configure the device ready for this experiment.

set_IA_precision();
% if (enable_default && any(strcmpi(p.UsingDefaults, 'IA_precision'))) || any(strcmpi(varargin, 'IA_precision'))
%     ziDAQ('setInt', ['/' device '/system/impedance/precision'], p.Results.IA_precision);
%     switch ziDAQ('getInt', ['/' device '/system/impedance/precision'])
%     case 0
%         IA_precision = 'low->fast';
%     case 1
%         IA_precision = 'high->medium';
%     case 2
%         IA_precision = 'very high->slow';
%     end
%     if major.disp, fprintf('IA precision set to %s.\n', IA_precision); end
%     SettingsStr.IA_precision = IA_precision;
% end
set_voltage_range();
% if (enable_default && any(strcmpi(p.UsingDefaults, 'voltage_range'))) || any(strcmpi(varargin, 'voltage_range'))
%     ziDAQ('setDouble', ['/' device '/sigouts/' c.out_c '/range'], p.Results.voltage_range);
%     if major.disp, fprintf('Signal out voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/sigouts/' c.out_c '/range'])); end
% 
%     ziDAQ('setDouble', ['/' device '/sigins/' c.in_c '/range'], p.Results.voltage_range);
%     if major.disp, fprintf('Signal in voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/sigins/' c.in_c '/range'])); end
% end

if (enable_default && any(strcmpi(p.UsingDefaults, 'imp50ohm'))) || any(strcmpi(varargin, 'imp50ohm'))
    ziDAQ('setInt', ['/' device '/sigins/' c.in_c '/imp50'], p.Results.imp50ohm)
    if ziDAQ('getInt', ['/' device '/sigins/' c.in_c '/imp50'])
    impohm = '50';
    else
    impohm = '10M';
    end
    if major.disp, fprintf('Output impedance set to %s ohm.\n', impohm); end
    SettingsStr.impohm = impohm;
end

ziDAQ('setInt', ['/' device '/sigouts/' c.out_c '/on'], 1);

ziDAQ('setDouble', ['/' device '/sigouts/' c.out_c '/enables/' c.out_mixer_c], 1); 

ziDAQ('setDouble', ['/' device '/demods/*/phaseshift'], 0);
if major.disp, fprintf('Demod phase shift set to %g.\n', ziDAQ('getDouble', ['/' device '/demods/' c.demod_c '/phaseshift'])); end

if (enable_default && any(strcmpi(p.UsingDefaults, 'demod_rate'))) || any(strcmpi(varargin, 'demod_rate'))
    ziDAQ('setDouble', ['/' device '/demods/' c.demod_c '/rate'], p.Results.demod_rate);
    if major.disp, fprintf('Demod data transfer rate set to %g Sa/sec. \n', ziDAQ('getDouble', ['/' device '/demods/' c.demod_c '/rate'])); end
end

ziDAQ('setInt', ['/' device '/demods/' c.demod_c '/harmonic'], 1);
if major.disp, fprintf('Demod harmonic set to %g.\n', ziDAQ('getInt', ['/' device '/demods/' c.demod_c '/harmonic'])); end

ziDAQ('setInt', ['/' device '/demods/' c.demod_c '/enable'], 1);
ziDAQ('setInt', ['/' device '/demods/*/oscselect'], str2double(c.osc_c));
ziDAQ('setInt', ['/' device '/demods/*/adcselect'], str2double(c.in_c));

if (enable_default && any(strcmpi(p.UsingDefaults, 'demod_time_constant'))) || any(strcmpi(varargin, 'demod_time_constant'))
    ziDAQ('setDouble', ['/' device '/demods/*/timeconstant'], p.Results.demod_time_constant);
    if major.disp, fprintf('Demod time constant set to %g sec.\n', ziDAQ('getDouble', ['/' device '/demods/' c.demod_c '/timeconstant'])); end
end

if (enable_default && any(strcmpi(p.UsingDefaults, 'two_terminal'))) || any(strcmpi(varargin, 'two_terminal'))
    ziDAQ('setInt', ['/' device '/imps/' c.imp_c '/mode'], p.Results.two_terminal);
    if ziDAQ('getInt', ['/' device '/imps/' c.imp_c '/mode'])
        Terminals = 'Two (2) Terminal';
    else
        Terminals = 'Four (4) Terminal';
    end
    if major.disp, fprintf('Measurement mode set to %s.\n', Terminals); end
    SettingsStr.Terminals = Terminals;
end


if ziDAQ('getInt', ['/' device '/imps/' c.imp_c '/mode'])
    if (enable_default && any(strcmpi(p.UsingDefaults, 'cable_length'))) || any(strcmpi(varargin, 'cable_length'))
        ziDAQ('setInt', ['/' device '/system/impedance/calib/cablelength'], p.Results.cable_length);
        if major.disp, fprintf('Cable length compensation set to %g m.\n', ziDAQ('getInt', ['/' device '/system/impedance/calib/cablelength'])); end
    end
    if (enable_default && any(strcmpi(p.UsingDefaults, 'AC'))) || any(strcmpi(varargin, 'AC'))
        ziDAQ('setInt', ['/' device '/sigins/' c.in_c '/ac'], p.Results.AC);
        if ziDAQ('getInt', ['/' device '/sigins/' c.in_c '/ac'])
            AC2T = 'ENABLED';
        else
            AC2T = 'DISABLED';
        end
        if major.disp, fprintf('High-pass filter for two terminal set to %s.\n', AC2T); end
        SettingsStr.AC = AC2T;
    end
else
    if (enable_default && any(strcmpi(p.UsingDefaults, 'AC'))) || any(strcmpi(varargin, 'AC'))
        ziDAQ('setInt', ['/' device '/imps/' c.imp_c '/ac'], p.Results.AC);
        if ziDAQ('getInt', ['/' device '/imps/' c.imp_c '/ac'])
            AC4T = 'ENABLED';
        else
            AC4T = 'DISABLED';
        end
        if major.disp, fprintf('High-pass filter for four terminal set to %s.\n', AC4T); end
        SettingsStr.AC = AC4T;
    end
end
ziDAQ('setInt', ['/' device '/imps/' c.imp_c '/demod/order'], p.Results.demod_LFP_order);
if major.disp, fprintf('IA Demod LFP order set to %g.\n', ziDAQ('getInt', ['/' device '/imps/' c.imp_c '/demod/order'])); end

if (enable_default && any(strcmpi(p.UsingDefaults, 'auto_range'))) || any(strcmpi(varargin, 'auto_range'))
    ziDAQ('setInt', ['/' device '/imps/' c.imp_c '/auto/inputrange'], p.Results.auto_range);
    if ziDAQ('getInt', ['/' device '/imps/' c.imp_c '/auto/inputrange'])
        auto_r = 'ENABLED';
    else
        auto_r = 'DISABLED';
    end
    if major.disp, fprintf('IA auto range set to %s.\n', auto_r); end
    SettingsStr.auto_r = auto_r;
end

set_current_range();
% if (enable_default && any(strcmpi(p.UsingDefaults, 'current_range'))) || any(strcmpi(varargin, 'current_range'))
%     ziDAQ('setDouble', ['/' device '/imps/' c.imp_c '/current/range'], p.Results.current_range);
%     if major.disp, fprintf('IA current range set to %g A.\n', ziDAQ('getDouble', ['/' device '/imps/' c.imp_c '/current/range'])); end
% end

% if (enable_default && any(strcmpi(p.UsingDefaults, 'voltage_range'))) || any(strcmpi(varargin, 'voltage_range'))
%     ziDAQ('setDouble', ['/' device '/imps/' c.imp_c '/voltage/range'], p.Results.voltage_range);
%     if major.disp, fprintf('IA input voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/imps/' c.imp_c '/voltage/range'])); end
% 
%     ziDAQ('setDouble', ['/' device '/imps/' c.imp_c '/output/range'], p.Results.voltage_range);
%     if major.disp, fprintf('IA output voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/imps/' c.imp_c '/output/range'])); end
% end

ziDAQ('setInt', ['/' device '/imps/' c.imp_c '/auto/output'], 0);
ziDAQ('setInt', ['/' device '/imps/' c.imp_c '/enable'], 1);
ziDAQ('setInt', ['/' device '/imps/' c.imp_c '/output/on'], 1);
    % 1 - Rs Cs, 2 - Rs Ls, 3 - G B, 4 -D Cs, 5 - Qs Cs, 6 - D Ls,
    % 7 - Q Ls, 8 - Rp Lp, 9 - D Cp
if (enable_default && any(strcmpi(p.UsingDefaults, 'model'))) || any(strcmpi(varargin, 'model'))
    ziDAQ('setInt', ['/' device '/imps/' c.imp_c '/model'], p.Results.model);
    switch ziDAQ('getInt', ['/' device '/imps/' c.imp_c '/model'])
        case 0
            model = 'Rp Cp';
        case 1
            model = 'Rs Cs';
        case 2
            model = 'Rs Ls';
        case 3
            model = 'G B';
        case 4
            model = 'D Cs';
        case 5
            model = 'Qs Cs';
        case 6
            model = 'D Ls';
        case 7
            model = 'Q Ls';
        case 8
            model = 'Rp Lp';
        case 9
            model = 'D Cp';
    end
    if major.disp, fprintf('IA parameter representation set to %s.\n', model); end
    SettingsStr.model = model;
end
OutSettings.str = SettingsStr;
actual_frequency_vec = zeros(1,length(frequency_vec));
actual_amplitude_vec = zeros(1,length(amplitude_vec));
actual_offset_vec = zeros(1,length(offset_vec));
title_params = '';
offset_sign = sign(offset_vec(round(length(offset_vec)/4)));
%% Sweep by
if additional_settings_internal.display.graph.save.if
    StartTime = datestr(now, 'yyyy-mm-dd hh-MM');
    SavePath = [additional_settings_internal.display.graph.save.path '\sweeps ' StartTime];
    mkdir(SavePath)
end
for v = 1:lpairs
    frequency_val = frequency_vec(min(v,lf));
    amplitude_val = amplitude_vec(min(v,la));
    offset_val = offset_vec(min(v,lo));
    if additional_settings_internal.output.DLCP
        offset_val = offset_val + offset_sign*amplitude_val;
    end
    if any(strcmpi(sweep_order(2:3), 'frequency'))
        ziDAQ('setDouble', ['/' device '/imps/' c.imp_c '/freq'], frequency_val)
        actual_frequency = ziDAQ('getDouble', ['/' device '/imps/' c.imp_c '/freq']);
        actual_frequency_vec(v) = actual_frequency;
        if major.disp, fprintf('IA frequency set to %g Hz.\n', actual_frequency_vec(v)); end
        title_params = [' Freq = ' num2str(actual_frequency)];
        for ic = size(SetByRange,2)
            if strcmpi(SetByRange{2,ic}, 'frequency')
                if (frequency_val>=SetByRange{3,ic}(1) && frequency_val<=SetByRange{3,ic}(2)) && SetByRange{5,ic}
                    eval(['set_' SetByRange{1,ic} '(' num2str(SetByRange{4,ic}) ');'])
                    SetByRange{5,ic} = false;
                elseif (frequency_val<SetByRange{3,ic}(1) || frequency_val>SetByRange{3,ic}(2)) && ~SetByRange{5,ic}
                    eval(['set_' SetByRange{1,ic} '();'])
                    SetByRange{5,ic} = true;
                end
            end
        end        
    end
    if any(strcmpi(sweep_order(2:3), 'amplitude'))
        ziDAQ('setInt', ['/' device '/imps/' c.imp_c '/auto/output'], 0);
        ziDAQ('setDouble', ['/' device '/imps/' c.imp_c '/output/amplitude'], amplitude_val)
        actual_amplitude = ziDAQ('getDouble', ['/' device '/imps/' c.imp_c '/output/amplitude']);
        actual_amplitude_vec(v) = actual_amplitude;
        if major.disp, fprintf('IA test signal (amplitude) set to %g mV.\n', 1000*actual_amplitude_vec(v)); end
        title_params = [title_params ' amp = ' num2str(actual_amplitude)];
        for ic = size(SetByRange,2)
            if strcmpi(SetByRange{2,ic}, 'amplitude') && SetByRange{5,ic}
                if (amplitude_val>=SetByRange{3,ic}(1) && amplitude_val<=SetByRange{3,ic}(2)) && SetByRange{5,ic}
                    eval(['set_' SetByRange{1,ic} '(' num2str(SetByRange{4,ic}) ');'])
                    SetByRange{5,ic} = false;
                elseif (amplitude_val<SetByRange{3,ic}(1) || amplitude_val>SetByRange{3,ic}(2)) && ~SetByRange{5,ic}
                    eval(['set_' SetByRange{1,ic} '();'])
                    SetByRange{5,ic} = true;
                end
            end
        end 
    end
    if any(strcmpi(sweep_order(2:3), 'offset'))
        ziDAQ('setDouble', ['/' device '/imps/' c.imp_c '/bias/value'], offset_val)
        ziDAQ('setInt', ['/' device '/imps/' c.imp_c '/bias/enable'], 1);
        actual_offset = ziDAQ('getDouble', ['/' device '/imps/' c.imp_c '/bias/value']);
        actual_offset_vec(v) = actual_offset;
        if major.disp, fprintf('IA bias voltage (offset) set to %g V.\n', actual_offset_vec(v)); end
        title_params = [title_params ' bias = ' num2str(actual_offset)];
        for ic = size(SetByRange,2)
            if strcmpi(SetByRange{2,ic}, 'offset')
                if (offset_val>=SetByRange{3,ic}(1) && offset_val<=SetByRange{3,ic}(2)) && SetByRange{5,ic}
                    eval(['set_' SetByRange{1,ic} '(' num2str(SetByRange{4,ic}) ');'])
                    SetByRange{5,ic} = false;
                elseif (offset_val<SetByRange{3,ic}(1) || offset_val>SetByRange{3,ic}(2)) && ~SetByRange{5,ic}
                    eval(['set_' SetByRange{1,ic} '();'])
                    SetByRange{5,ic} = true;
                end
            end
        end 
    end
    %% Perform Parameter (=sweep_order{1}) Sweep and save fig if enabled
    [select_data_one, full_data_one, sw_plot] = MFIA_general_sweeper(device, additional_settings_internal, sweep_order(1), sweep_range, pts, read_param_struct, unmatched_vars{:});
    if additional_settings_internal.display.graph.save.if
        title(title_params);
        saveas(sw_plot,[SavePath '\' title_params '.png']);
        title_params = '';
    end
    %%
    % Disable repetitive text display
    if major.disp && ~additional_settings_internal.display.text.major.each_sweep
        additional_settings_internal.display.text.major.disp = false;
        major.disp = false;
        fprintf('Major text display will not be shown henceforth.\n')
    end
    if minor.disp && ~additional_settings_internal.display.text.minor.each_sweep
        additional_settings_internal.display.text.minor.disp = false;
        minor.disp = false;
        fprintf('Minor text display (settings) will not be shown henceforth.\n')
    end
    
    for stc = sweep_order(2:3)
        st = stc{:};
        eval(['select_data_one.actualVals.' st '=actual_' st '_vec(v);'])
        eval(['full_data_one.actualVals.' st '=actual_' st '_vec(v);'])
        eval(['select_data_one.inputVals.' st '=' st '_vec(v);'])
        eval(['full_data_one.inputVals.' st '=' st '_vec(v);'])
    end

    
    select_data(v) = select_data_one;
    full_data(v) = full_data_one;
end

% Disable everything when the run is finished
ziDisableEverything(device);
%% Function
function Get = set_IA_precision(varargin)
    if (enable_default && any(strcmpi(p.UsingDefaults, 'IA_precision'))) || any(strcmpi(vararginTOP, 'IA_precision'))
        if isempty(varargin)
            SetVar = p.Results.IA_precision;
        else
            SetVar = varargin{1};
        end
        ziDAQ('setInt', ['/' device '/system/impedance/precision'], SetVar);
        Get = ziDAQ('getInt', ['/' device '/system/impedance/precision']);
        switch Get 
        case 0
            IA_precision = 'low->fast';
        case 1
            IA_precision = 'high->medium';
        case 2
            IA_precision = 'very high->slow';
        end
        if major.disp, fprintf('IA precision set to %s.\n', IA_precision); end
        if isempty(varargin)
            SettingsStr.IA_precision = IA_precision;
        end
    end
end

function Get = set_voltage_range(varargin)
    if (enable_default && any(strcmpi(p.UsingDefaults, 'voltage_range'))) || any(strcmpi(vararginTOP, 'voltage_range'))
        if isempty(varargin)
            SetVar = p.Results.voltage_range;
        else
            SetVar = varargin{1};
        end
        ziDAQ('setDouble', ['/' device '/sigouts/' c.out_c '/range'], SetVar);
        if major.disp, fprintf('Signal out voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/sigouts/' c.out_c '/range'])); end

        ziDAQ('setDouble', ['/' device '/sigins/' c.in_c '/range'], SetVar);
        if major.disp, fprintf('Signal in voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/sigins/' c.in_c '/range'])); end
        
        ziDAQ('setDouble', ['/' device '/imps/' c.imp_c '/voltage/range'], SetVar);
        Get = ziDAQ('getDouble', ['/' device '/imps/' c.imp_c '/voltage/range']);
        if major.disp, fprintf('IA input voltage range set to %g V.\n', Get); end

        ziDAQ('setDouble', ['/' device '/imps/' c.imp_c '/output/range'], SetVar);
        if major.disp, fprintf('IA output voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/imps/' c.imp_c '/output/range'])); end
    end
end
function Get = set_current_range(varargin)
    if (enable_default && any(strcmpi(p.UsingDefaults, 'current_range'))) || any(strcmpi(vararginTOP, 'current_range'))
        if isempty(varargin)
            SetVar = p.Results.current_range;
        else
            SetVar = varargin{1};
        end
        ziDAQ('setDouble', ['/' device '/imps/' c.imp_c '/current/range'], SetVar);
        Get = ziDAQ('getDouble', ['/' device '/imps/' c.imp_c '/current/range']);
        if major.disp, fprintf('IA current range set to %g A.\n', Get); end
    end
end
end
