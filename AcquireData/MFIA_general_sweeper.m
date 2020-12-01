function [select_data, full_data, varargout]  = MFIA_general_sweeper(device, additional_settings, sweep_param, sweep_range, pts, read_param_struct, intermediate_read, varargin)
varargout = {};
% Define parameters relevant to this example. Default values specified by the
% inputParser below are overwritten if specified as name-value pairs via the
% `varargin` input argument.
p = inputParser;
p.KeepUnmatched=true;
isnonnegscalar = @(x) isnumeric(x) && isscalar(x) && (x > 0);

%Sweep timeout.
p.addParameter('timeout', 120, isnonnegscalar);

%Fetch data during the sweep.
p.addParameter('intermediate_read', 0, @isnumeric);

% Perform one single sweep.
p.addParameter('loopcount', 1, isnonnegscalar);
% Logarithmic sweep mode.
p.addParameter('xmapping', 0, @isnumeric);
% Binary scan type.
p.addParameter('scan', 1, @isnumeric);


% The value used for the Sweeper's 'settling/inaccuracy' parameter: This
% defines the settling time the sweeper should wait before changing a sweep
% parameter and recording the next sweep data point. The settling time is
% calculated from the specified proportion of a step response function that
% should remain. The value provided here, 0.001, is appropriate for fast and
% reasonably accurate amplitude measurements. For precise noise measurements
% it should be set to ~100n.
% Note: The actual time the sweeper waits before
% recording data is the maximum time specified by settling/time and
% defined by settling/inaccuracy.
p.addParameter('sweep_inaccuracy', 0.001, @isnumeric);
% We don't require a fixed settling/time since there is no DUT involved
% in this example's setup (only a simple feedback cable) so set this to
% zero. We need only wait for the filter response to settle, specified via
% settling/inaccuracy.
p.addParameter('settling_time', 0, @isnumeric);
% Minimum time to record and average data is 50 time constants.
p.addParameter('averaging_time_constant', 50, @isnumeric);
% Minimal number of samples that we want to record and average is 100. Note,
% the number of samples used for averaging will be the maximum number of
% samples specified by either averaging/tc or averaging/sample.
p.addParameter('averaging_samples', 100, isnonnegscalar);



% Use automatic bandwidth control for each measurement.
% For fixed bandwidth, set bandwidthcontrol to 1 and specify a bandwidth.
% For manual bandwidth control, set  bandwidthcontrol to 2. bandwidth must also be set
% to a value > 0 although it is ignored. Otherwise Auto control is automatically chosen (for backwards compatibility reasons).
% ziDAQ('set', h, 'bandwidth', 100);
p.addParameter('bandwidth_control', 2, @isnumeric);
% Sets the bandwidth overlap mode (default 0). If ENABLED, the bandwidth of a
% sweep point may overlap with the frequency of neighboring sweep points. The
% effective bandwidth is only limited by the maximal bandwidth setting and
% omega suppression. As a result, the bandwidth is independent of the number
% of sweep points. For frequency response analysis bandwidth overlap should be
% ENABLED to achieve maximal sweep speed (default: 0). 0 = Disable, 1 = Enable.
p.addParameter('bandwidth_overlap', 0, @isnumeric);

p.parse(varargin{:});


keys = {'demodulator', 'demod', 'impedance', 'imp'};
values = {'demod', 'demod', 'impedance', 'impedance'};
node_dic = containers.Map(keys,values);
read_param_cell = fn_struct2cell(read_param_struct);
grid_search = strfind(read_param_cell(1,:), 'grid');
if ~any([grid_search{:}])
    read_param_cell = [{'read_param_struct.demod.grid'; true} read_param_cell];
end
if ~isempty(read_param_cell)
    for i = 1:size(read_param_cell,2)
        for st = {'demod', 'imp'}
            locs = strfind(read_param_cell{1,i},'.');
            if contains(read_param_cell{1,i}(locs(1)+1:locs(2)-1),st{:})
            read_param_cell{1,i} = ['select_data' read_param_cell{3,i}];
            read_param_cell{5,i} = ['sample_' node_dic(st{:}) read_param_cell{1,i}(locs(2):end)];
            read_param_cell{6,i} = regexprep(read_param_cell{1,i}(8:end), {'\.' '(' ')'},{'_', '' ''});
            end
        end
    end
end
           
       
%% Set default additional settings
% Define device channels.
additional_settings_internal.channels.demod_c = '0'; % demod channel, for paths on the device
additional_settings_internal.channels.demod_idx = str2double(additional_settings_internal.channels.demod_c)+1; % 1-based indexing, to access the data
additional_settings_internal.channels.out_c = '0'; % signal output channel
% Get the value of the instrument's default Signal Output mixer channel.
additional_settings_internal.channels.out_mixer_c = '1';
additional_settings_internal.channels.in_c = '0'; % signal input channel
additional_settings_internal.channels.osc_c = '0'; % oscillator
additional_settings_internal.channels.imp_c = '0'; % IA channel
additional_settings_internal.channels.imp_index = str2double(additional_settings_internal.channels.imp_c)+1; % IA, 1-based indexing, to access the data
% Graphs
additional_settings_internal.display.graph.disp = true;
additional_settings_internal.display.graph.during_sweep = false;
% Text
additional_settings_internal.display.text.major.disp = true;
additional_settings_internal.display.text.minor.disp = true;

% Overwrite default additional settings
additional_settings_internal = update_structure(additional_settings_internal, additional_settings); 

channels = additional_settings_internal.channels;
channels_cell = fn_struct2cell(channels);

for i = 1:size(channels_cell,2)
    eval([channels_cell{4,i} '=' channels_cell{1,i} ';']);
end


major.disp = additional_settings_internal.display.text.major.disp;
minor.disp = additional_settings_internal.display.text.minor.disp;
graph.disp = additional_settings_internal.display.graph.disp;
graph.during_sweep = additional_settings_internal.display.graph.during_sweep;

%% Sweeper settings
% Create a thread for the sweeper
h = ziDAQ('sweep');
% Device on which sweeping will be performed
ziDAQ('set', h, 'device', device);
if major.disp, fprintf('Device set to %s.\n', device); end
%% Set sweep parameters

% Perform sweeps consisting of sweep_samplecount measurement points (i.e.,
% record the subscribed data for sweep_samplecount different frequencies).
ziDAQ('set', h, 'samplecount', pts);

if strcmp(sweep_param, 'frequency')
    % Sweeping setting is the frequency of the output signal
    ziDAQ('set', h, 'gridnode', ['oscs/' osc_c '/freq']); 
    % Start frequency = 1 kHz
    ziDAQ('set', h, 'start', sweep_range(1));
    % Stop frequency
    ziDAQ('set', h, 'stop', sweep_range(2));
    % Plot function
    plot_func = @semilogx;
    % Plot label
    lbl = 'Oscillator Frequency [Hz]';
	if major.disp, fprintf('Frequency Sweep from %g Hz to %g Hz, %d pts.\n', ziDAQ('get', h, 'start').start, ziDAQ('get', h, 'stop').stop, ziDAQ('get', h, 'samplecount').samplecount); end
elseif strcmp(sweep_param, 'amplitude')
    % Sweeping setting is the amplitude of the output signal
    ziDAQ('set', h, 'gridnode', ['sigouts/' out_c '/amplitudes/' out_mixer_c]);
    % Start test signal
    ziDAQ('set', h, 'start', sweep_range(1));
    % Stop test signal
    ziDAQ('set', h, 'stop', sweep_range(2));
    % Plot function
    plot_func = @plot;
    % Plot label
    lbl = 'Test Signal Amplitude [V]';
	if major.disp, fprintf('Test Signal Sweep from %.2f mV to %.2f mV, %d pts.\n', 1000*ziDAQ('get', h, 'start').start, 1000*ziDAQ('get', h, 'stop').stop, ziDAQ('get', h, 'samplecount').samplecount); end
elseif strcmp(sweep_param, 'offset')
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
end

ziDAQ('set', h, 'loopcount', p.Results.loopcount);
if minor.disp, fprintf('Loop count set to %d.\n', ziDAQ('get', h, 'loopcount').loopcount); end

ziDAQ('set', h, 'xmapping', p.Results.xmapping);
if ziDAQ('get', h, 'xmapping').xmapping
    xmapping = 'Linear';
else
    xmapping = 'Logarithmic';
end
if minor.disp, fprintf('X axis mapping set to %s.\n', xmapping); end

ziDAQ('set', h, 'scan', p.Results.scan);
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

ziDAQ('set', h, 'averaging/tc', p.Results.averaging_time_constant);
if minor.disp, fprintf('Minimum time to record and average data set to %g time constants.\n', ziDAQ('get', h, 'averaging/tc').averaging.tc); end

ziDAQ('set', h, 'averaging/sample', p.Results.averaging_samples);
if minor.disp, fprintf('Averaging set to %g samples.\n', ziDAQ('get', h, 'averaging/sample').averaging.sample); end

ziDAQ('set', h, 'bandwidthcontrol', p.Results.bandwidth_control);
switch ziDAQ('get', h, 'bandwidthcontrol').bandwidthcontrol
    case 0 
        bandwidth_control = 'Manual';
    case 1 
        bandwidth_control = 'Fixed';
	case 2 
        bandwidth_control = 'Auto';
end
if minor.disp, fprintf('Bandwidth control set to %s.\n', bandwidth_control); end

ziDAQ('set', h, 'bandwidthoverlap', p.Results.bandwidth_overlap);
if ziDAQ('get', h, 'bandwidthoverlap').bandwidthoverlap
    bandwidth_overlap = 'DISABLED';
else
    bandwidth_overlap = 'ENABLED';
end
if minor.disp, fprintf('Bandwidth overlap set to %s.\n', bandwidth_overlap); end


% Subscribe to the node from which data will be recorded.
ziDAQ('subscribe', h, ['/' device '/demods/' demod_c '/sample']);
ziDAQ('subscribe', h, ['/' device, '/imps/' imp_c '/sample']);

% Start sweeping.
ziDAQ('execute', h);
if major.disp, fprintf('Sweep Started\n'); end

full_data = [];
select_data = struct;

figure(1); clf;
t0 = tic;
% Read and plot intermediate data until the sweep has finished.
while ~ziDAQ('finished', h)
    pause(1);
    if intermediate_read
        tmp = ziDAQ('read', h);
        if major.disp, fprintf('Sweep progress %0.0f%%\n', ziDAQ('progress', h) * 100); end
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
                        eval(c{1})
                        eval([varargout{i} '=' c{5} ';'])
                        i=i+1;
                    end
                end
                % Plot data during sweep
                if graph.disp && graph.during_sweep
                % Sweep parameter values at which measurement points were taken
                    sweep_param_arr = sample_impedance.grid;
        %             valid = ~isnan(sweep_i);
                     plot_data(plot_func, lbl, sweep_param_arr, select_data, read_param_cell{1,:})
        %             drawnow;
                end
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
                    eval(c{1})
                    eval([varargout{i} '=' c{5} ';'])
                    i=i+1;
                end
            end
        end
        % Plot data
        if graph.disp
        % Frequency values at which measurement points were taken
        sweep_param_arr = sample_impedance.grid;
        % Plot the final result
        plot_data(plot_func, lbl, sweep_param_arr, select_data, read_param_cell{1,:})
        end
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

function plot_data(plot_func, lbl, x, select_data, varargin)
% Plot data
clf
[subrow, subcol] = subplot_min_rectangle(length(varargin));
i=1;
for c = varargin
    subplot(subrow, subcol, i)
    s = plot_func(x, eval(c{:}));
    set(s, 'LineWidth', 1.5)
    set(s, 'Color', 'black');
    grid on
    xlabel(lbl)
    ylabel(c{:})
    i = i+1;
end
end
