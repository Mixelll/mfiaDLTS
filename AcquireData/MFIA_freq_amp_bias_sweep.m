function [select_data, full_data, varargout]  = MFIA_freq_amp_bias_sweep(device_id, additional_settings, sweep_order, sweep_range, pts, frequency_vec, amplitude_vec, offset_vec, read_param_struct, varargin)
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
lf = length(frequency_vec);
la = length(amplitude_vec);
lo = length(offset_vec);
if ~(sum([lf la lo]>0)>=2 && (lf==la || lf==lo || lo==la))
    fprintf('sweep_order #2 and #3 vecs do not match in length');
    return
end
    
clear ziDAQ;

if ~exist('device_id', 'var')
    error(['No value for device_id specified. The first argument to the ' ...
           'example should be the device ID on which to run the example, ' ...
           'e.g. ''dev2006'' or ''uhf-dev2006''.'])
end

% Check the ziDAQ MEX (DLL) and Utility functions can be found in Matlab's path.
if ~(exist('ziDAQ') == 3) && ~(exist('ziCreateAPISession', 'file') == 2)
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
if ~any(strcmpi([branches], 'DEMODS'))
  fprintf('\nThis example requires lock-in functionality which is not available on %s.\n', device);
  return
end

% Define parameters relevant to this example. Default values specified by the
% inputParser below are overwritten if specified as name-value pairs via the
% `varargin` input argument.
p = inputParser;
p.KeepUnmatched=true;
isnumcalar = @(x) isnumeric(x) && isscalar(x) && (x > 0);
isnonnegscalar = @(x) isnumeric(x) && isscalar(x) && (x > 0);

% Set IA parameter extraction model (from impedance). 0 - Rp Cp,
% 1 - Rs Cs, 2 - Rs Ls, 3 - G B, 4 -D Cs, 5 - Qs Cs, 6 - D Ls,
% 7 - Q Ls, 8 - Rp Lp, 9 - D Cp
p.addParameter('model', 0, @isnumcalar);

% Enable two terminal measurement.
p.addParameter('two_terminal', 1, @isnumcalar);

% Enable two terminal cable length compensation [m].
p.addParameter('cable_length', 1, @isnumcalar);

% Enable high-pass filter.
p.addParameter('AC', 0, @isnumcalar);

% Enable four terminal high-pass filter.
p.addParameter('imp50ohm', 0, @issnumcalar);

% Enable current auto range.
p.addParameter('auto_range', 0, isnumcalar);

% Current manual range, [A].
p.addParameter('current_range', 10e-6, @isnumeric);

% Voltage manual range, [V].
p.addParameter('voltage_range', 3, @isnumeric);

% demod time constant, [s].
p.addParameter('demod_time_constant', 0.007, @isnumeric);

% demod rate, [Hz].
p.addParameter('demod_rate', 13e3, @isnumeric);

p.parse(varargin{:});
UsingDefaults = p.UsingDefaults;
unmatched_vars = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)];
unmatched_vars = unmatched_vars.';

%% Set default additional settings
additional_settings_internal.enable_default = true;
% Define device channels.
additional_settings_internal.channels.demod_c = '0'; % demod channel, for paths on the device
additional_settings_internal.channels.demod_idx = str2double(additional_settings_internal.channels.demod_c)+1; % 1-based indexing, to access the data
additional_settings_internal.channels.out_c = '0'; % signal output channel
% Get the value of the instrument's default Signal Output mixer channel.
additional_settings_internal.channels.out_mixer_c = num2str(ziGetDefaultSigoutMixerChannel(props, str2num(additional_settings_internal.channels.out_c)));
additional_settings_internal.channels.in_c = '0'; % signal input channel
additional_settings_internal.channels.osc_c = '0'; % oscillator
additional_settings_internal.channels.imp_c = '0'; % IA channel
additional_settings_internal.channels.imp_index = str2double(additional_settings_internal.channels.imp_c)+1; % IA, 1-based indexing, to access the data
% Graphs
additional_settings_internal.display.graph.disp = true;
additional_settings_internal.display.graph.during_sweep = false;
% Text
additional_settings_internal.display.text.major.disp = true;
additional_settings_internal.display.text.major.each_sweep = true;
additional_settings_internal.display.text.minor.disp = true;
additional_settings_internal.display.text.minor.each_sweep = false;

% Overwrite default additional settings
additional_settings_internal = update_structure(additional_settings_internal, additional_settings); 

channels = additional_settings_internal.channels;
channels_cell = fn_struct2cell(channels);

for i = 1:size(channels_cell,2)
    eval([channels_cell{4,i} '=' channels_cell{1,i} ';']);
end

major.disp = additional_settings_internal.display.text.major.disp;
enable_default = additional_settings_internal.enable_default;

% Create a base configuration: Disable all available outputs, awgs,
% demods, scopes,...
ziDisableEverything(device);

%% Configure the device ready for this experiment.

if enable_default && any(strcmpi(UsingDefaults, 'voltage_range'))
    ziDAQ('setDouble', ['/' device '/sigouts/' out_c '/range'], p.Results.voltage_range);
    if major.disp, fprintf('Signal out voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/sigouts/' out_c '/range'])); end

    ziDAQ('setDouble', ['/' device '/sigins/' in_c '/range'], p.Results.voltage_range);
    if major.disp, fprintf('Signal in voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/sigins/' in_c '/range'])); end
end

if enable_default && any(strcmpi(UsingDefaults, 'imp50ohm'))
    ziDAQ('setInt', ['/' device '/sigins/' in_c '/imp50'], p.Results.imp50ohm)
    if ziDAQ('getInt', ['/' device '/sigins/' in_c '/imp50'])
    impohm = '50';
    else
    impohm = '10M';
    end
    if major.disp, fprintf('Output impedance set to %s ohm.\n', impohm); end
end

ziDAQ('setInt', ['/' device '/sigouts/' out_c '/on'], 1);

ziDAQ('setDouble', ['/' device '/sigouts/' out_c '/enables/' out_mixer_c], 1); 

ziDAQ('setDouble', ['/' device '/demods/*/phaseshift'], 0);
if major.disp, fprintf('Demod phase shift set to %g.\n', ziDAQ('getDouble', ['/' device '/demods/' demod_c '/phaseshift'])); end

ziDAQ('setInt', ['/' device '/demods/*/order'], 4);
if major.disp, fprintf('Demod order set to %g.\n', ziDAQ('getInt', ['/' device '/demods/' demod_c '/order'])); end

if enable_default && any(strcmpi(UsingDefaults, 'demod_rate'))
    ziDAQ('setDouble', ['/' device '/demods/' demod_c '/rate'], p.Results.demod_rate);
    if major.disp, fprintf('Demod rate set to %g. Hz\n', ziDAQ('getDouble', ['/' device '/demods/' demod_c '/rate'])); end
end

ziDAQ('setInt', ['/' device '/demods/' demod_c '/harmonic'], 1);
if major.disp, fprintf('Demod harmonic set to %g.\n', ziDAQ('getInt', ['/' device '/demods/' demod_c '/harmonic'])); end

ziDAQ('setInt', ['/' device '/demods/' demod_c '/enable'], 1);
ziDAQ('setInt', ['/' device '/demods/*/oscselect'], str2double(osc_c));
ziDAQ('setInt', ['/' device '/demods/*/adcselect'], str2double(in_c));

if enable_default && any(strcmpi(UsingDefaults, 'demod_time_constant'))
    ziDAQ('setDouble', ['/' device '/demods/*/timeconstant'], p.Results.demod_time_constant);
    if major.disp, fprintf('Demod time constant set to %g sec.\n', ziDAQ('getDouble', ['/' device '/demods/' demod_c '/timeconstant'])); end
end

if enable_default && any(strcmpi(UsingDefaults, 'two_terminal'))
    ziDAQ('setInt', ['/' device '/imps/' imp_c '/mode'], p.Results.two_terminal);
    if ziDAQ('getInt', ['/' device '/imps/' imp_c '/mode'])
        Tselect = 'Two (2) Terminal';
    else
        Tselect = 'Four (4) Terminal';
    end
    if major.disp, fprintf('Measurement mode set to %s.\n', Tselect); end
end

if ziDAQ('getInt', ['/' device '/imps/' imp_c '/mode'])
    if enable_default && any(strcmpi(UsingDefaults, 'cable_length'))
        ziDAQ('setInt', ['/' device '/system/impedance/calib/cablelength'], p.Results.cable_length);
        if major.disp, fprintf('Cable length compensation set to %g m.\n', ziDAQ('getInt', ['/' device '/system/impedance/calib/cablelength'])); end
    end
    if enable_default && any(strcmpi(UsingDefaults, 'AC'))
        ziDAQ('setInt', ['/' device '/sigins/' in_c '/ac'], p.Results.AC);
        if ziDAQ('getInt', ['/' device '/sigins/' in_c '/ac'])
            AC2T = 'ENABLED';
        else
            AC2T = 'DISABLED';
        end
        if major.disp, fprintf('High-pass filter for two terminal set to %s.\n', AC2T); end
    end
else
    if enable_default && any(strcmpi(UsingDefaults, 'AC'))
        ziDAQ('setInt', ['/' device '/imps/' imp_c '/ac'], p.Results.AC);
        if ziDAQ('getInt', ['/' device '/imps/' imp_c '/ac'])
            AC4T = 'ENABLED';
        else
            AC4T = 'DISABLED';
        end
        if major.disp, fprintf('High-pass filter for four terminal set to %s.\n', AC4T); end
    end
end

if enable_default && any(strcmpi(UsingDefaults, 'auto_range'))
    ziDAQ('setInt', ['/' device '/imps/' imp_c '/auto/inputrange'], p.Results.auto_range);
    if ziDAQ('getInt', ['/' device '/imps/' imp_c '/mode'])
        auto_r = 'ENABLED';
    else
        auto_r = 'DISABLED';
    end
    if major.disp, fprintf('IA auto range set to %s.\n', auto_r); end
end

if enable_default && any(strcmpi(UsingDefaults, 'current_range'))
    ziDAQ('setDouble', ['/' device '/imps/' imp_c '/current/range'], p.Results.current_range);
    if major.disp, fprintf('IA current range set to %g A.\n', ziDAQ('getDouble', ['/' device '/imps/' imp_c '/current/range'])); end
end

if enable_default && any(strcmpi(UsingDefaults, 'voltage_range'))
    ziDAQ('setDouble', ['/' device '/imps/' imp_c '/voltage/range'], p.Results.voltage_range);
    if major.disp, fprintf('IA input voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/imps/' imp_c '/voltage/range'])); end

    ziDAQ('setDouble', ['/' device '/imps/' imp_c '/output/range'], p.Results.voltage_range);
    if major.disp, fprintf('IA output voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/imps/' imp_c '/output/range'])); end
end

ziDAQ('setInt', ['/' device '/imps/' imp_c '/auto/output'], 0);
ziDAQ('setInt', ['/' device '/imps/' imp_c '/enable'], 1);
ziDAQ('setInt', ['/' device '/imps/' imp_c '/output/on'], 1);

if enable_default && any(strcmpi(UsingDefaults, 'model'))
    ziDAQ('setInt', ['/' device '/imps/' imp_c '/model'], p.Results.model);
    if major.disp, fprintf('IA output voltage range set to %g V.\n', ziDAQ('getInt', ['/' device '/imps/' imp_c '/model'])); end
end

actual_frequency_vec = zeros(1,length(frequency_vec));
actual_amplitude_vec = zeros(1,length(amplitude_vec));
actual_offset_vec = zeros(1,length(offset_vec));

%% Sweep by
for v = 1:max([lf,la,lo])
    if any(strcmpi(sweep_order(2:3), 'frequency'))
        ziDAQ('setDouble', ['/' device '/imps/' imp_c '/freq'], frequency_vec(v))
        actual_frequency_vec(v) = ziDAQ('getDouble', ['/' device '/imps/' imp_c '/freq']);
        if major.disp, fprintf('IA frequency set to %g Hz.\n', actual_frequency_vec(v)); end
    end
    if any(strcmpi(sweep_order(2:3), 'amplitude'))
        ziDAQ('setInt', ['/' device '/imps/' imp_c '/auto/output'], 0);
        ziDAQ('setDouble', ['/' device '/imps/' imp_c '/output/amplitude'], amplitude_vec(v))
        actual_amplitude_vec(v) = ziDAQ('getDouble', ['/' device '/imps/' imp_c '/output/amplitude']);
        if major.disp, fprintf('IA test signal (amplitude) set to %g mV.\n', 1000*actual_amplitude_vec(v)); end
    end
    if any(strcmpi(sweep_order(2:3), 'offset'))
        ziDAQ('setDouble', ['/' device '/imps/' imp_c '/bias/value'], offset_vec(v))
        ziDAQ('setInt', ['/' device '/imps/' imp_c '/bias/enable'], 1);
        actual_offset_vec(v) = ziDAQ('getDouble', ['/' device '/imps/' imp_c '/bias/value']);
        if major.disp, fprintf('IA bias voltage (offset) set to %g V.\n', actual_offset_vec(v)); end
    end
    
    [select_data_one, full_data_one] = MFIA_general_sweeper(device, additional_settings_internal, sweep_order(1), sweep_range, pts, read_param_struct, unmatched_vars{:});
    if ~additional_settings_internal.display.text.major.each_sweep
        additional_settings_internal.display.text.major.disp = false;
        major.disp = false;
        fprintf('Major text display will not be shown henceforth.\n')
    end
    if ~additional_settings_internal.display.text.minor.each_sweep
        additional_settings_internal.display.text.minor.disp = false;
        fprintf('Minor text display (settings) will not be shown henceforth.\n')
    end
    
    for st = sweep_order(2:3)
        st = st{:};
        eval(['select_data_one.' st '=actual_' st '_vec(v);'])
        eval(['full_data_one.' st '=actual_' st '_vec(v);'])
    end
    
    select_data(v) = select_data_one;
    full_data(v) = full_data_one;
end


% Sweeper module returns a structure with following elements:
% * timestamp -> Time stamp data [uint64]. Divide the timestamp by the
% device's clockbase in order to get seconds, the clockbase can be obtained
% via: clockbase = double(ziDAQ('getInt', ['/' device '/clockbase']));
% * x -> Demodulator x value in volt [double]
% * y -> Demodulator y value in volt [double]
% * r -> Demodulator r value in Vrms [double]
% * phase ->  Demodulator theta value in rad [double]
% * xstddev -> Standard deviation of demodulator x value [double]
% * ystddev -> Standard deviation of demodulator x value [double]
% * rstddev -> Standard deviation of demodulator r value [double]
% * phasestddev -> Standard deviation of demodulator theta value [double]
% * grid ->  Values of sweeping setting (frequency values at which
% demodulator samples where recorded) [double]
% * bandwidth ->  Filter bandwidth for each measurement point [double].
% * tc ->  Filter time constant for each measurement point [double].
% * settling ->  Waiting time for each measurement point [double]
% * frequency ->  Oscillator frequency for each measurement point [double]
% (in his case same as grid).
% * frequencystddev -> Standard deviation of oscillator frequency
% * realbandwidth -> The actual bandwidth set
% Assuming we're doing a frequency sweep:
% * settimestamp -> The time at which we verify that the frequency for the
% current sweep point was set on the device (by reading back demodulator data)
% * nexttimestamp -> The time at which we can get the data for that sweep point.
% i.e., nexttimestamp - settimestamp corresponds roughly to the
% demodulator filter settling time.

end
