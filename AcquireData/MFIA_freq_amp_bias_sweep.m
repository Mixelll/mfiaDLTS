function [select_data, full_data, varargout]  = MFIA_freq_amp_bias_sweep(device_id, device_properties, sweep_order, sweep_range, pts, frequency_vec, amplitude_vec, offset_vec, read_param_struct, intermediate_read, varargin)
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
if ~any(strcmp([branches], 'DEMODS'))
  fprintf('\nThis example requires lock-in functionality which is not available on %s.\n', device);
  return
end

% Define parameters relevant to this example. Default values specified by the
% inputParser below are overwritten if specified as name-value pairs via the
% `varargin` input argument.
p = inputParser;
p.KeepUnmatched=true;
isnonnegscalar = @(x) isnumeric(x) && isscalar(x) && (x > 0);

% Enable two terminal measurement.
p.addParameter('two_terminal', 1, @isnumeric);

% Enable two terminal cable length compensation [m].
p.addParameter('cable_length', 1, @isnumeric);

% Enable four terminal high-pass filter.
p.addParameter('AC', 0, @isnumeric);

% Enable current auto range.
p.addParameter('auto_range', 0, @isnumeric);

% Current manual range, [A].
p.addParameter('current_range', 10e-6, @isnumeric);

% Voltage manual range, [V].
p.addParameter('voltage_range', 3, @isnumeric);

% demod time constant, [s].
p.addParameter('demod_time_constant', 0.007, @isnumeric);

% demod rate, [Hz].
p.addParameter('demod_rate', 13e3, @isnumeric);

p.parse(varargin{:});

unmatched_vars = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)];
unmatched_vars = unmatched_vars.';

% Define some other helper parameters.
demod_c = '0'; % demod channel, for paths on the device
demod_idx = str2double(demod_c)+1; % 1-based indexing, to access the data
out_c = '0'; % signal output channel
% Get the value of the instrument's default Signal Output mixer channel.
out_mixer_c = num2str(ziGetDefaultSigoutMixerChannel(props, str2num(out_c)));
in_c = '0'; % signal input channel
osc_c = '0'; % oscillator
imp_c = '0'; % IA channel
imp_index = str2double(imp_c)+1; % IA, 1-based indexing, to access the data

if isa(device_properties, 'struct') && any(strcmp(fieldnames(device_properties), 'channels'))
    for c = fieldnames(device_properties.channels)
        eval([c{:} '=device_properties.channels.' c{:}]);
    end
end

% Create a base configuration: Disable all available outputs, awgs,
% demods, scopes,...
ziDisableEverything(device);

%% Configure the device ready for this experiment.

ziDAQ('setDouble', ['/' device '/sigouts/' out_c '/range'], p.Results.voltage_range);
ziDAQ('setInt', ['/' device '/sigins/' in_c '/imp50'], 1);
ziDAQ('setDouble', ['/' device '/sigins/' in_c '/range'], p.Results.voltage_range);
ziDAQ('setInt', ['/' device '/sigouts/' out_c '/on'], 1);

% ziDAQ('setDouble', ['/' device '/sigouts/' out_c '/amplitudes/*'], 0);
% ziDAQ('setDouble', ['/' device '/sigouts/' out_c '/amplitudes/' out_mixer_c], p.Results.amplitude);
ziDAQ('setDouble', ['/' device '/sigouts/' out_c '/enables/' out_mixer_c], 1); 
ziDAQ('setDouble', ['/' device '/demods/*/phaseshift'], 0);
ziDAQ('setInt', ['/' device '/demods/*/order'], 4);
ziDAQ('setDouble', ['/' device '/demods/' demod_c '/rate'], p.Results.demod_rate);
ziDAQ('setInt', ['/' device '/demods/' demod_c '/harmonic'], 1);
ziDAQ('setInt', ['/' device '/demods/' demod_c '/enable'], 1);
ziDAQ('setInt', ['/' device '/demods/*/oscselect'], str2double(osc_c));
ziDAQ('setInt', ['/' device '/demods/*/adcselect'], str2double(in_c));
ziDAQ('setDouble', ['/' device '/demods/*/timeconstant'], p.Results.demod_time_constant);
% ziDAQ('setDouble', ['/' device '/oscs/' osc_c '/freq'], ); % [Hz]

ziDAQ('setInt', ['/' device '/imps/' imp_c '/mode'], p.Results.two_terminal);
if p.Results.two_terminal
    ziDAQ('setInt', ['/' device '/system/impedance/calib/cablelength'], p.Results.cable_length);
else
    ziDAQ('setInt', ['/' device '/imps/' imp_c '/ac'], p.Results.AC);
    ziDAQ('setInt', ['/' device '/sigins/' in_c '/ac'], p.Results.AC);
end

ziDAQ('setInt', ['/' device '/imps/' imp_c '/auto/inputrange'], p.Results.auto_range);
ziDAQ('setDouble', ['/' device '/imps/' imp_c '/current/range'], p.Results.current_range);
ziDAQ('setDouble', ['/' device '/imps/' imp_c '/voltage/range'], p.Results.voltage_range);
ziDAQ('setDouble', ['/' device '/imps/' imp_c '/output/range'], p.Results.voltage_range);
ziDAQ('setInt', ['/' device '/imps/' imp_c '/auto/output'], 0);
ziDAQ('setInt', ['/' device '/imps/' imp_c '/enable'], 1);
ziDAQ('setInt', ['/' device '/imps/' imp_c '/output/on'], 1);


%% Sweep by
select_data = repmat(struct,1,max([lf la lo]));
full_data = select_data;
figure(1); clf;
for v = 1:max([lf,la,lo])
    if any(strcmp(sweep_order(2:3), 'frequency'))
        ziDAQ('setDouble', ['/' device 'imps/' imp_c '/bias/value'], frequency_vec(v)) 
    end
    if any(strcmp(sweep_order(2:3), 'amplitude'))
        ziDAQ('setInt', ['/' device '/imps/' imp_c '/auto/output'], 0);
        ziDAQ('setDouble', ['/' device 'imps/' imp_c '/outputs/amplitude'], amplitude_vec(v))
    end
    if any(strcmp(sweep_order(2:3), 'offset'))
        ziDAQ('setDouble', ['/' device 'imps/' imp_c '/bias/value'], offset_vec(v))
        ziDAQ('setInt', ['/' device '/imps/' imp_c '/bias/enable'], 1);
    end
    
    [select_data_one, full_data_one] = MFIA_general_sweeper(device, device_properties, sweep_order(1), sweep_range, pts, read_param_struct, intermediate_read, unmatched_vars{:});
    
    for st = sweep_order(2:3)
        st = st{:};
        eval(['select_data_one.' st '=' st '_vec(v);'])
        eval(['full_data_one.' st '=' st '_vec(v);'])
    end
    
    select_data(v) = select_data_one;
    full_data(v) = full_data_one;
end

% Release module resources. Especially important if modules are created
% inside a loop to prevent excessive resource consumption.
ziDAQ('clear', h);

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

function plot_data(frequencies, r, theta, amplitude, style)
% Plot data
clf
subplot(2, 1, 1)
s = semilogx(frequencies, 20*log10(r*2*sqrt(2)/amplitude), style);
set(s, 'LineWidth', 1.5)
set(s, 'Color', 'black');
grid on
xlabel('Frequency [Hz]')
ylabel('Amplitude [dBV]')
subplot(2, 1, 2)
s = semilogx(frequencies, theta*180/pi, style);
set(s, 'LineWidth', 1.5)
set(s, 'Color', 'black');
grid on
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
end
