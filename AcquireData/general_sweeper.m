function data = general_sweeper(device, sweep_param, sweep_range, pts, intermediate_read, timeout)

% Define parameters relevant to this example. Default values specified by the
% inputParser below are overwritten if specified as name-value pairs via the
% `varargin` input argument.
p = inputParser;
isnonnegscalar = @(x) isnumeric(x) && isscalar(x) && (x > 0);

% The value used for the Sweeper's 'samplecount' parameter: This
% specifies the number of points that will be swept (i.e., the number of
% frequencies swept in a frequency sweep).
p.addParamValue('length', 100, isnonnegscalar);

% Perform one single sweep.
p.addParamValue('loopcount', 1, @isnumeric);
% Logarithmic sweep mode.
p.addParamValue('xmapping', 0, @isnumeric);
% Binary scan type.
p.addParamValue('scan', 1, @isnumeric);



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
p.addParamValue('sweep_inaccuracy', 0.001, @isnumeric);
% We don't require a fixed settling/time since there is no DUT involved
% in this example's setup (only a simple feedback cable) so set this to
% zero. We need only wait for the filter response to settle, specified via
% settling/inaccuracy.
p.addParamValue('settling_time', 0, @isnumeric);
% Minimum time to record and average data is 50 time constants.
p.addParamValue('averaging_time_constant', 50, @isnumeric);
% Minimal number of samples that we want to record and average is 100. Note,
% the number of samples used for averaging will be the maximum number of
% samples specified by either averaging/tc or averaging/sample.
p.addParamValue('averaging_samples', 100, @isnumeric);



% Use automatic bandwidth control for each measurement.
% For fixed bandwidth, set bandwidthcontrol to 1 and specify a bandwidth.
% For manual bandwidth control, set  bandwidthcontrol to 2. bandwidth must also be set
% to a value > 0 although it is ignored. Otherwise Auto control is automatically chosen (for backwards compatibility reasons).
% ziDAQ('set', h, 'bandwidth', 100);
p.addParamValue('bandwidth_control', 2, @isnumeric);
% Sets the bandwidth overlap mode (default 0). If enabled, the bandwidth of a
% sweep point may overlap with the frequency of neighboring sweep points. The
% effective bandwidth is only limited by the maximal bandwidth setting and
% omega suppression. As a result, the bandwidth is independent of the number
% of sweep points. For frequency response analysis bandwidth overlap should be
% enabled to achieve maximal sweep speed (default: 0). 0 = Disable, 1 = Enable.
p.addParamValue('bandwidth_overlap', 0.1, @isnumeric);

p.parse(varargin{:});


% Define some other helper parameters.
demod_c = '0'; % demod channel, for paths on the device
demod_idx = str2double(demod_c)+1; % 1-based indexing, to access the data
out_c = '0'; % signal output channel
% Get the value of the instrument's default Signal Output mixer channel.
out_mixer_c = num2str(ziGetDefaultSigoutMixerChannel(props, str2num(out_c)));
in_c = '0'; % signal input channel
osc_c = '0'; % oscillator

if isa(device_properties, 'struct') && any(strcmp(fieldnames(device_properties), 'chennels'))
    for c = fieldnames(device_properties.channels)
        eval([c{:} '=device_properties.channels.' c{:}]);
    end
end




%% Sweeper settings
% Create a thread for the sweeper
h = ziDAQ('sweep');
% Device on which sweeping will be performed
ziDAQ('set', h, 'device', device);
%% Sweep by
if strcamp(sweep_param, 'osc_frequency')
    % Sweeping setting is the frequency of the output signal
    ziDAQ('set', h, 'gridnode', ['oscs/' osc_c '/freq']); 
    % Start frequency = 1 kHz
    ziDAQ('set', h, 'start', sweep_range(1));
    % Stop frequency
    ziDAQ('set', h, 'stop', sweep_range(2));
	fprintf('Frequency Sweep from %g Hz to %g Hz, %d pts.\n', sweep_range(1), sweep_range(2), pts);
elseif strcamp(sweep_param, 'test_signal')
    % Sweeping setting is the amplitude of the output signal
    ziDAQ('set', h, 'gridnode', ['sigouts/' out_c '/amplitudes/' out_mixer_c]);
    % Start test signal
    ziDAQ('set', h, 'start', sweep_range(1));
    % Stop test signal
    ziDAQ('set', h, 'stop', sweep_range(2));
	fprintf('Test Signal Sweep from %.2f mV to %.2f mV, %d pts.\n', 1000*sweep_range(1), 1000*sweep_range(2), pts);
elseif strcamp(sweep_param, 'bias_voltage')
    % Sweeping setting is the offset of the output signal
    ziDAQ('set', h, 'gridnode', ['sigouts/' osc_c '/offset']);
    % Start bias voltage
    ziDAQ('set', h, 'start', sweep_range(1));
    % Stop bias voltage
    ziDAQ('set', h, 'stop', sweep_range(2));
	fprintf('Bias Voltage Sweep from %.2f V to %.2f V, %d pts.\n', sweep_range(1), sweep_range(2), pts);
end

% Perform sweeps consisting of sweep_samplecount measurement points (i.e.,
% record the subscribed data for sweep_samplecount different frequencies).
ziDAQ('set', h, 'samplecount', pts);

ziDAQ('set', h, 'loopcount', p.Results.loopcount);
fprintf('Loop count set to %d  pts.\n', p.Results.loopcount);

ziDAQ('set', h, 'xmapping', p.Results.xmapping);
if p.Results.xmapping
    xmapping = 'Linear';
else
    xmapping = 'Logarithmic';
end
fprintf('X axis mapping set to %s.\n', xmapping);

ziDAQ('set', h, 'scan', p.Results.scan);
switch p.Results.scan
    case 0 
        scan = 'Sequential';
    case 1 
        scan = 'Binary';
	case 2 
        scan = 'Bidirectional';
	case 3 
        scan = 'Reverse';
end
fprintf('Scan type set to %s.\n', scan);

ziDAQ('set', h, 'settling/time', p.Results.settling_time);
fprintf('settling/time set to %g sec.\n', p.Results.settling_time);

ziDAQ('set', h, 'settling/inaccuracy', p.Results.sweep_inaccuracy);
fprintf('Sweep inaccuracy set to %g.\n', p.Results.sweep_inaccuracy);

ziDAQ('set', h, 'averaging/tc', p.Results.averaging_time_constant);
fprintf('Minimum time to record and average data set to %g time constants.\n', p.Results.averaging_time_constant);

ziDAQ('set', h, 'averaging/sample', p.Results.averaging_samples);
fprintf('Averaging set to %g samples.\n', p.Results.averaging_time_constant);

ziDAQ('set', h, 'bandwidthcontrol', p.Results.bandwidth_control);
switch p.Results.bandwidth_control
    case 0 
        bandwidth_control = 'Manual';
    case 1 
        bandwidth_control = 'Fixed';
	case 2 
        bandwidth_control = 'Auto';
end
fprintf('Bandwidth control set to %s.\n', bandwidth_control);

if p.Results.xmapping
    bandwidth_overlap = 'Disabled';
else
    bandwidth_overlap = 'Enabled';
end
ziDAQ('set', h, 'bandwidthoverlap', p.Results.bandwidth_overlap);
fprintf('Bandwidth overlap set to %s.\n', bandwidth_overlap);


% Subscribe to the node from which data will be recorded.
ziDAQ('subscribe', h, ['/' device '/demods/' demod_c '/sample']);

% Start sweeping.
ziDAQ('execute', h);
fprintf('Sweep Started\n');

data = [];
eval([sweep_param '_i= nan(1, pts)']);
r = nan(1, pts);
theta = nan(1, pts);
p0 = nan(1, pts);
p1 = nan(1, pts);

figure(1); clf;
t0 = tic;
% Read and plot intermediate data until the sweep has finished.
while ~ziDAQ('finished', h)
    pause(1);
    if intermediate_read
        tmp = ziDAQ('read', h);
        fprintf('Sweep progress %0.0f%%\n', ziDAQ('progress', h) * 100);
        % Using intermediate reads we can plot a continuous refinement of the ongoing
        % measurement. If not required it can be removed.
        if ziCheckPathInData(tmp, ['/' device '/demods/' demod_c '/sample'])
            sample = tmp.(device).demods(demod_idx).sample{1};
            if ~isempty(sample)
                data = tmp;
                % Get the magnitude and phase of demodulator from the sweeper result.
                r = sample.r;
                theta = sample.phase;
                param0 = sample.param0;
                param1 = sample.param1;
                % Frequency values at which measurement points were taken
                eval([sweep_param '_i= sample.grid']);
    %             valid = ~isnan(sweep_i);
    %             plot_data(sweep_i(valid), r(valid), theta(valid), p.Results.amplitude, '.-')
    %             drawnow;
            end
        end
    end
    if toc(t0) > timeout
        ziDAQ('clear', h);
        error('Timeout: Sweeper failed to finish after %f seconds.', timeout)
    end
end
fprintf('Sweep completed after %.2f s.\n', toc(t0));


% Read the data. This command can also be executed during the waiting (as above).
tmp = ziDAQ('read', h);

% Unsubscribe from the node; stop filling the data from that node to the
% internal buffer in the Data Server.
ziDAQ('unsubscribe', h, ['/' device '/demods/*/sample']);

% Process any remainging data returned by read().
if ziCheckPathInData(tmp, ['/' device '/demods/' demod_c '/sample'])
    sample = tmp.(device).demods(demod_idx).sample{1};
    if ~isempty(sample)
        % Extracting R component and phase of input signal
        % As several sweeps may be returned, a cell array is used.
        % In this case we pick the first sweep result by {1}.
        data = tmp;
        r = sample.r;
        theta = sample.phase;
        param0 = sample.param0;
        param1 = sample.param1;
        % Frequency values at which measurement points were taken
        eval([sweep_param '_i= sample.grid']);
        % Plot the final result
%         plot_data(sweep_i, r, theta, p.Results.amplitude, '-')
        
    end
end


% Release module resources. Especially important if modules are created
% inside a loop to prevent excessive resource consumption.
ziDAQ('clear', h);

for s = {'r', 'theta', 'param0', 'param1', [sweep_param '_i']}
    eval(['data.' s '=' s])
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
