function [IV_data, select_data, full_data, sw_plot, varargout]  = MFIA_DC(device_id, additional_settings, sweep_range, pts, read_param_struct, InSavePath, hysteresis, varargin)
varargout = {};

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

% Enable four terminal high-pass filter.
p.addParameter('imp50ohm', 0, isnumscalar);

% Enable current auto range.
p.addParameter('auto_range', 1, isnumscalar);

% Current manual range, [A].
p.addParameter('current_range', 10e-6, @isnumeric);

% Voltage manual range, [V].
p.addParameter('voltage_range', 3, @isnumeric);

% demod time constant, [s].
p.addParameter('demod_time_constant', 0.007, @isnumeric);

% demod data transfer rate, [Sa/sec].
p.addParameter('demod_rate', 13.39e3, @isnumeric);


%Sweep timeout.
p.addParameter('timeout', 120, isnonnegscalar);
% Perform one single sweep.
p.addParameter('loopcount', 1, isnonnegscalar);
% Logarithmic sweep mode.

% Binary scan type.
p.addParameter('scan', 1, @isnumeric);

% Minimum wait time in seconds between a sweep parameter change and the recording of the next sweep point. This
% parameter can be used to define the required settling time of the experimental setup. The effective wait time
% is the maximum of this value and the demodulator filter settling time determined from the Inaccuracy value specified
p.addParameter('settling_time', 0, @isnumeric);
% Demodulator filter settling inaccuracy defining the wait time between a sweep parameter change and
% recording of the next sweep point. The settling time is calculated as the time required to attain the specified
% remaining proportion [1e-13,0.1] of an incoming step function. Typical inaccuracy the number of filter time
% constants the sweeper has to wait. The maximum between this value and the settling time is taken as wait time until the
% next sweep point is recorded values: 10 m for highest sweep speed for large signals, 100 u for precise amplitude
% measurements, 100 n for precise noise measurements. Depending on the order the settling accuracy will define
p.addParameter('sweep_inaccuracy', 0.001, @isnumeric);
% Sets the number of data samples per sweeper parameter point that is considered in the measurement. The maximum
% between samples, time and number of time constants is taken as effective calculation time.
p.addParameter('averaging_samples', 20, isnonnegscalar);
% Sets the time during which data samples are processed. The maximum between samples, time and number
% of time constants is taken as effective calculation time
p.addParameter('averaging_time', 0.1, @isnumeric);
% Sets the effective measurement time per sweeper parameter point that is considered in the
% measurement. The maximum between samples, time and number of time constants is
% taken as effective calculation time.
p.addParameter('averaging_time_constant', 15, @isnumeric);

p.parse(varargin{:});

keys = {'demodulator', 'demod', 'impedance', 'imp'};
values = {'demod', 'demod', 'impedance', 'impedance'};
node_dic = containers.Map(keys,values);

read_param_cell = fn_struct2cell(read_param_struct);
grid_search = strfind(read_param_cell(1,:), 'grid');
if ~any([grid_search{:}])
    read_param_cell = [{'read_param_struct.demod.grid'; true; '.demod.grid'; grid} read_param_cell];
end
for i = 1:size(read_param_cell,2)
    read_param_cell{1,i} = ['select_data' read_param_cell{3,i}];
    for st = {'demod', 'imp'}
        locs = strfind(read_param_cell{1,i},'.');
        if contains(read_param_cell{1,i}(locs(1)+1:locs(2)-1),st{:})
        read_param_cell{5,i} = ['sample_' node_dic(st{:}) read_param_cell{1,i}(locs(2):end)];
        read_param_cell{6,i} = regexprep(read_param_cell{1,i}(8:end), {'\.' '(' ')'},{'_', '' ''});
        end
    end
end
           
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

channels = additional_settings_internal.channels;
channels_cell = fn_struct2cell(channels);

for i = 1:size(channels_cell,2)
    eval([channels_cell{4,i} '=' channels_cell{1,i} ';']);
end

major.disp = additional_settings_internal.display.text.major.disp;
minor.disp = additional_settings_internal.display.text.minor.disp;
enable_default = additional_settings_internal.enable_default;

% Create a base configuration: Disable all available outputs, awgs, demods, scopes,...
ziDisableEverything(device);       

major.disp = additional_settings_internal.display.text.major.disp;
minor.disp = additional_settings_internal.display.text.minor.disp;
graph.disp = additional_settings_internal.display.graph.disp;
graph.during_sweep = additional_settings_internal.display.graph.during_sweep;


if ~isempty(InSavePath)
    StartTime = datestr(now, 'yyyy-mm-dd hh-MM');
    SavePath0 = [InSavePath '\DC'];
    i=1;
    while true
        SavePath = [SavePath0 num2str(i)];
        if ~exist(SavePath, 'dir')
            mkdir(SavePath)
            break
        else
            i=i+1;
        end
    end
end



%% Configure the device ready for this experiment.

ziDAQ('setDouble', ['/' device '/imps/' imp_c '/freq'], 0)
freq = ziDAQ('getDouble', ['/' device '/imps/' imp_c '/freq']);
if major.disp, fprintf('IA frequency set to %g Hz.\n', freq); end


if (enable_default && any(strcmpi(p.UsingDefaults, 'voltage_range'))) || any(strcmpi(varargin, 'voltage_range'))
    ziDAQ('setInt', ['/' device '/system/impedance/precision'], p.Results.IA_precision);
    switch ziDAQ('getInt', ['/' device '/system/impedance/precision'])
    case 0
        IA_precision = 'low->fast';
    case 1
        IA_precision = 'high->medium';
    case 2
        IA_precision = 'very high->slow';
    end
    if major.disp, fprintf('IA precision set to %s.\n', IA_precision); end
end

if (enable_default && any(strcmpi(p.UsingDefaults, 'voltage_range'))) || any(strcmpi(varargin, 'voltage_range'))
    ziDAQ('setDouble', ['/' device '/sigouts/' out_c '/range'], p.Results.voltage_range);
    if major.disp, fprintf('Signal out voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/sigouts/' out_c '/range'])); end

    ziDAQ('setDouble', ['/' device '/sigins/' in_c '/range'], p.Results.voltage_range);
    if major.disp, fprintf('Signal in voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/sigins/' in_c '/range'])); end
end

if (enable_default && any(strcmpi(p.UsingDefaults, 'imp50ohm'))) || any(strcmpi(varargin, 'imp50ohm'))
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

if (enable_default && any(strcmpi(p.UsingDefaults, 'demod_rate'))) || any(strcmpi(varargin, 'demod_rate'))
    ziDAQ('setDouble', ['/' device '/demods/' demod_c '/rate'], p.Results.demod_rate);
    if major.disp, fprintf('Demod data transfer rate set to %g Sa/sec. \n', ziDAQ('getDouble', ['/' device '/demods/' demod_c '/rate'])); end
end

ziDAQ('setInt', ['/' device '/demods/' demod_c '/enable'], 1);
ziDAQ('setInt', ['/' device '/demods/*/oscselect'], str2double(osc_c));
ziDAQ('setInt', ['/' device '/demods/*/adcselect'], str2double(in_c));

if (enable_default && any(strcmpi(p.UsingDefaults, 'demod_time_constant'))) || any(strcmpi(varargin, 'demod_time_constant'))
    ziDAQ('setDouble', ['/' device '/demods/*/timeconstant'], p.Results.demod_time_constant);
    if major.disp, fprintf('Demod time constant set to %g sec.\n', ziDAQ('getDouble', ['/' device '/demods/' demod_c '/timeconstant'])); end
end

if (enable_default && any(strcmpi(p.UsingDefaults, 'two_terminal'))) || any(strcmpi(varargin, 'two_terminal'))
    ziDAQ('setInt', ['/' device '/imps/' imp_c '/mode'], p.Results.two_terminal);
    if ziDAQ('getInt', ['/' device '/imps/' imp_c '/mode'])
        Tselect = 'Two (2) Terminal';
    else
        Tselect = 'Four (4) Terminal';
    end
    if major.disp, fprintf('Measurement mode set to %s.\n', Tselect); end
end

if (enable_default && any(strcmpi(p.UsingDefaults, 'auto_range'))) || any(strcmpi(varargin, 'auto_range'))
    ziDAQ('setInt', ['/' device '/imps/' imp_c '/auto/inputrange'], p.Results.auto_range);
    if ziDAQ('getInt', ['/' device '/imps/' imp_c '/mode'])
        auto_r = 'ENABLED';
    else
        auto_r = 'DISABLED';
    end
    if major.disp, fprintf('IA auto range set to %s.\n', auto_r); end
end

if (enable_default && any(strcmpi(p.UsingDefaults, 'current_range'))) || any(strcmpi(varargin, 'current_range'))
    ziDAQ('setDouble', ['/' device '/imps/' imp_c '/current/range'], p.Results.current_range);
    if major.disp, fprintf('IA current range set to %g A.\n', ziDAQ('getDouble', ['/' device '/imps/' imp_c '/current/range'])); end
end

if (enable_default && any(strcmpi(p.UsingDefaults, 'voltage_range'))) || any(strcmpi(varargin, 'voltage_range'))
    ziDAQ('setDouble', ['/' device '/imps/' imp_c '/voltage/range'], p.Results.voltage_range);
    if major.disp, fprintf('IA input voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/imps/' imp_c '/voltage/range'])); end

    ziDAQ('setDouble', ['/' device '/imps/' imp_c '/output/range'], p.Results.voltage_range);
    if major.disp, fprintf('IA output voltage range set to %g V.\n', ziDAQ('getDouble', ['/' device '/imps/' imp_c '/output/range'])); end
end

ziDAQ('setInt', ['/' device '/imps/' imp_c '/auto/output'], 0);
ziDAQ('setInt', ['/' device '/imps/' imp_c '/enable'], 1);
ziDAQ('setInt', ['/' device '/imps/' imp_c '/output/on'], 1);
    % 1 - Rs Cs, 2 - Rs Ls, 3 - G B, 4 -D Cs, 5 - Qs Cs, 6 - D Ls,
    % 7 - Q Ls, 8 - Rp Lp, 9 - D Cp
if (enable_default && any(strcmpi(p.UsingDefaults, 'model'))) || any(strcmpi(varargin, 'model'))
    ziDAQ('setInt', ['/' device '/imps/' imp_c '/model'], p.Results.model);
    switch ziDAQ('getInt', ['/' device '/imps/' imp_c '/model'])
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
end


%% Sweeper settings
% Create a thread for the sweeper
h = ziDAQ('sweep');
% Device on which sweeping will be performed
ziDAQ('set', h, 'device', device);
if minor.disp, fprintf('Device set to %s.\n', device); end

%% Set sweep parameters

% Perform sweeps consisting of sweep_samplecount measurement points (i.e.,
% record the subscribed data for sweep_samplecount different frequencies).
ziDAQ('set', h, 'samplecount', pts);

% Sweeping setting is the offset of the output signal
ziDAQ('set', h, 'gridnode', ['sigouts/' osc_c '/offset']);
% Start bias voltage
ziDAQ('set', h, 'start', sweep_range(1));
% Stop bias voltage
ziDAQ('set', h, 'stop', sweep_range(2));
% Plot function
plot_func = @plot;
% Plot label
lbl = 'Bias Voltage (offset) [V]';
if major.disp, fprintf('Bias Voltage Sweep from %.2f V to %.2f V, %d pts.\n', ziDAQ('get', h, 'start').start, ziDAQ('get', h, 'stop').stop, ziDAQ('get', h, 'samplecount').samplecount); end

ziDAQ('set', h, 'loopcount', p.Results.loopcount);
if minor.disp, fprintf('Loop count set to %d.\n', ziDAQ('get', h, 'loopcount').loopcount); end

ziDAQ('set', h, 'xmapping', 0);
if ziDAQ('get', h, 'xmapping').xmapping
    xmapping = 'Linear';
else
    xmapping = 'Logarithmic';
end
if minor.disp, fprintf('X axis mapping set to %s.\n', xmapping); end

if hysteresis
    scan = 2;
end
ziDAQ('set', h, 'scan', max(p.Results.scan, scan));
switch ziDAQ('get', h, 'scan').scan
    case 0 
        scan = 'Sequential';
    case 1 
        scan = 'Binary';
	case 2 
        scan = 'Bidirectional';
	case 3 
        scan = 'Reverse';
end
if minor.disp, fprintf('Scan type set to %s.\n', scan); end


ziDAQ('set', h, 'settling/time', p.Results.settling_time);
if minor.disp, fprintf('settling/time set to %g sec.\n', ziDAQ('get', h, 'settling/time').settling.time); end

ziDAQ('set', h, 'settling/inaccuracy', p.Results.sweep_inaccuracy);
if minor.disp, fprintf('Sweep inaccuracy set to %g.\n', ziDAQ('get', h, 'settling/inaccuracy').settling.inaccuracy); end

ziDAQ('set', h, 'averaging/sample', p.Results.averaging_samples);
if minor.disp, fprintf('Minimum averaging set to %g samples.\n', ziDAQ('get', h, 'averaging/sample').averaging.sample); end

ziDAQ('set', h, 'averaging/time', p.Results.averaging_time);
if minor.disp, fprintf('Minimum time to record and average data set to %g sec.\n', ziDAQ('get', h, 'averaging/time').averaging.time); end

ziDAQ('set', h, 'averaging/tc', p.Results.averaging_time_constant);
if minor.disp, fprintf('Minimum time to record and average data set to %g time constants.\n', ziDAQ('get', h, 'averaging/tc').averaging.tc); end


% Subscribe to the node from which data will be recorded.
ziDAQ('subscribe', h, ['/' device '/demods/' demod_c '/sample']);
ziDAQ('subscribe', h, ['/' device, '/imps/' imp_c '/sample']);

% Start sweeping.
ziDAQ('execute', h);
if major.disp, fprintf('Sweep Started\n'); end

full_data = [];
select_data = struct;
Title = [num2str(sweep_range(1)) ' to ' num2str(sweep_range(1))];

t0 = tic;
% Read and plot intermediate data until the sweep has finished.
while ~ziDAQ('finished', h)
    pause(1);
        if major.disp, fprintf('Sweep progress %0.0f%%\n', ziDAQ('progress', h) * 100); end
        % Plot data during sweep
        if graph.disp && graph.during_sweep
        tmp = ziDAQ('read', h);
        % Using intermediate reads we can plot a continuous refinement of the ongoing
        % measurement. If not required it can be removed.
            if ziCheckPathInData(tmp, ['/' device '/demods/' demod_c '/sample'])
                sample_demod = tmp.(device).demods(demod_idx).sample{1};
                sample_impedance = tmp.(device).imps(imp_index).sample{1};
                if ~isempty(sample_demod)
                    full_data = tmp;
                    % Get the desired parameters from the sweeper result.
                    i = 1;
                    for c = read_param_cell
                        if ~isempty(c{5}) && c{2}
                            eval([c{1} '=' c{5} ';'])
                            eval(['varargout{i} =' c{1} ';'])
                            i=i+1;
                        end
                    end
                    % Sweep parameter values at which measurement points were taken
                    sweep_param_arr = sample_impedance.grid;
                    sw_plot = plot_data(plot_func, lbl, sweep_param_arr, select_data, read_param_cell, Title, '');
%                     drawnow;
                end
            end
        end
    if toc(t0) > p.Results.timeout
        ziDAQ('clear', h);
        error('Timeout: Sweeper failed to finish after %f seconds.', p.Results.timeout)
    end
end
if major.disp, fprintf('Sweep completed after %.2f s.\n', toc(t0)); end


% Read the data. This command can also be executed during the waiting (as above).
tmp = ziDAQ('read', h);

% Unsubscribe from the node; stop filling the data from that node to the
% internal buffer in the Data Server.
ziDAQ('unsubscribe', h, ['/' device '/demods/*/sample']);
ziDAQ('unsubscribe', h, ['/' device, '/imps/' imp_c '/sample']);

% Process any remainging data returned by read().
    if ziCheckPathInData(tmp, ['/' device '/demods/' demod_c '/sample'])
        sample_demod = tmp.(device).demods(demod_idx).sample{1};
        sample_impedance = tmp.(device).imps(imp_index).sample{1};
        if ~isempty(sample_demod)
            full_data = tmp;
            % Get the desired parameters from the sweeper result.
            i = 1;
            for c = read_param_cell
                if ~isempty(c{5}) && c{2}
                    eval([c{1} '=' c{5} ';'])
                    eval(['varargout{i} =' c{1} ';'])
                    i=i+1;
                end
            end
        end
        % Plot data
        if graph.disp
        % Frequency values at which measurement points were taken
        sweep_param_arr = sample_impedance.grid;
        % Plot the final result
        sw_plot = plot_data(plot_func, lbl, sweep_param_arr, select_data, read_param_cell, Title, SavePath);
        end
    end
    
% Release module resources. Especially important if modules are created
% inside a loop to prevent excessive resource consumption.
ziDAQ('clear', h);
sweepfull.V = select_data.impedance.grid;
sweepfull.R = select_data.impedance.param0;
sweepfull.I = select_data.demod.r;
sweep.V = select_data.impedance.grid(1:pts);
sweep.R = select_data.impedance.param0(1:pts);
sweep.I = select_data.demod.r(1:pts);
if hysteresis
    sweeph.V = select_data.impedance.grid(pts+1:end);
    sweeph.R = select_data.impedance.param0(pts+1:end);
    sweeph.I = select_data.demod.r(pts+1:end);
end

IV_data.sweepfull = sweepfull;
IV_data.sweep = sweep;
if hysteresis
    IV_data.sweeph = sweeph;
end
save([SavePath '\IV'],'IV_data')
save([SavePath '\IV' StartTime],'IV_data') 
end

function s = plot_data(plot_func, lbl, x, select_data, read_param_cell, Title, savepath)
% Plot data
figure(1)
clf
subp_size = 0;
for c = read_param_cell
    if c{2} && ~contains(c{1}, 'grid')
        subp_size = subp_size+1;
    end
end
[subrow, subcol] = subplot_min_rectangle(subp_size);
i=1;
for c = read_param_cell
    if c{2} && ~contains(c{1}, 'grid')
        subplot(subrow, subcol, i)
        s = plot_func(x, eval(c{1}));
        set(s, 'LineWidth', 1.5)
        set(s, 'Color', 'black');
        grid on 
        xlabel(lbl)
        ylabel(c{4})
        i=i+1;
    end
end
title(Title);
if ~isempty(savepath)
    saveas(s,[savepath '\' Title '.png']);
    saveas(s,[savepath '\' Title '.fig']);
end
end
