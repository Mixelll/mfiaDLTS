%
% USAGE DATA = EXAMPLE_SWEEPER(DEVICE_ID)
%
% Perform a frequency/test_signal (amplitude)/bias_voltage (offset) sweep and gather demodulator data.
%
% Note that there are two name sets depending on context: 
% impedance analyzer context (sweeper context) - frequency (frequency), test_signal (amplitude), bias_voltage (offset)
%
% NOTE Please ensure that the ziDAQ folders 'Driver' and 'Utils' are in your
% Matlab path. To do this (temporarily) for one Matlab session please navigate
% to the ziDAQ base folder containing the 'Driver', 'Examples' and 'Utils'
% subfolders and run the Matlab function ziAddPath().
% >>> ziAddPath;
%
clear
% MFIA ID
device_id = 'dev5168';
desired_order = {'frequency', 'amplitude', 'offset'};
% Selected data to read (grid = sweep parameter)
% read_param_struct.demod.grid = true;
read_param_struct.demod.r = true;
read_param_struct.demod.phase = true;
read_param_struct.impedance.grid = true;
read_param_struct.impedance.param0 = true;
read_param_struct.impedance.param1 = true;

% Sweep by Frequency and iterate over AC amplitude and offset
sweep_order = {'frequency', 'amplitude', 'offset'};
start_frequency = 100; stop_frequency = 500e3; pts_frequency = 100; % Hz
start_amplitude = 0.05; stop_amplitude = 0.3; pts_amplitude = 2; % V
start_offset = 1; stop_offset = -1; pts_offset = 2; % V

[sweep_range, sweep_pts, frequency_vec, amplitude_vec, offset_vec] = MFIA_freq_amp_bias_value_pairs_withParser(sweep_order, 'start_frequency', start_frequency,...
    'stop_frequency', stop_frequency, 'pts_frequency', pts_frequency, 'start_amplitude', start_amplitude, 'stop_amplitude', stop_amplitude,...
    'pts_amplitude', pts_amplitude, 'start_offset', start_offset, 'stop_offset', stop_offset, 'pts_offset', pts_offset); 
[overwrite_defaults, additional_settings] = settings();

[select_data_sweep_order_struct_vec, full_data_sweep_order_struct_vec] = MFIA_freq_amp_bias_sweep(device_id, additional_settings, sweep_order, sweep_range, sweep_pts, frequency_vec, amplitude_vec, offset_vec, read_param_struct, overwrite_defaults{:});

[select_data_desired_order_3D, select_data_sweep_order_3D] = MFIA_data_reshape_3D(select_data_sweep_order_struct_vec, desired_order, sweep_order, pts_frequency, pts_amplitude, pts_offset, frequency_vec, amplitude_vec, offset_vec);


%% Overwite Defaults (uncomment and change value)
function [overwrite_defaults, additional_settings] = settings()
overwrite_defaults = {}; % don't touch
additional_settings = struct; % don't touch

    % BY SETTING THIS, THE INTERNAL SCRIPT DEFAULTS WILL OVERWRITE LABONE GUI INPUT
    % THE INTERNAL DEFAULTS ARE INSIDE MFIA_freq_amp_bias_sweep.
    % VALUES YOU UN-COMMENT BELOW WILL OVERWRITE REGARDLESS of "enable_default"
% additional_settings_internal.enable_default = true;
    % Graphs
% additional_settings.display.graph.disp = true;
% additional_settings.display.graph.during_sweep = true;
    % Text
% additional_settings.display.text.major.disp = true;
% additional_settings.display.text.major.each_sweep = true;
% additional_settings.display.text.minor.disp = true;
% additional_settings.display.text.minor.each_sweep = false;
    % Enable two terminal measurement.
% overwrite_defaults{:,end+1} = {'two_terminal'; 1};
    
    % Enable two terminal cable length compensation.
% overwrite_defaults{:,end+1} = {'cable_length'; 1};

    % Enable high-pass filter.
% overwrite_defaults{:,end+1} = {'AC'; 0};

    % Enable 50ohm output impedance. Disabled state is 10M ohm.
% overwrite_defaults{:,end+1} = {'imp50ohm'; 0};

    % Enable auto range.
% overwrite_defaults{:,end+1} = {'auto_range'; 0};

    % Current range.
% overwrite_defaults{:,end+1} = {'current_range'; 10e-6};

    % Voltage range.
% overwrite_defaults{:,end+1} = {'voltage_range'; 3};

    % Sweep timeout.
% overwrite_defaults{:,end+1} = {'timeout'; 120};

    % Fetch data during the sweep.
% overwrite_defaults{:,end+1} = {'intermediate_read'; 0};

    % Perform one single sweep.
% overwrite_defaults{:,end+1} = {'loopcount'; 1};

    % Logarithmic sweep mode.
% overwrite_defaults{:,end+1} = {'xmapping'; 0};

    % Binary scan type.
% overwrite_defaults{:,end+1} = {'scan'; 1};

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
%overwrite_defaults{:,end+1} = {'sweep_inaccuracy'; 0.001};

    % We don't require a fixed settling/time since there is no DUT involved
    % in this example's setup (only a simple feedback cable) so set this to
    % zero. We need only wait for the filter response to settle, specified via
    % settling/inaccuracy.
% overwrite_defaults{:,end+1} = {'settling_time'; 0};

    % Minimum time to record and average data is 50 time constants.	
% overwrite_defaults{:,end+1} = {'averaging_time_constant'; 50};

    % Minimal number of samples that we want to record and average is 100. Note,
    % the number of samples used for averaging will be the maximum number of
    % samples specified by either averaging/tc or averaging/sample.
% overwrite_defaults{:,end+1} = {'averaging_samples'; 100};

    % Use automatic bandwidth control for each measurement.
    % For fixed bandwidth, set bandwidthcontrol to 1 and specify a bandwidth.
    % For manual bandwidth control, set  bandwidthcontrol to 2. bandwidth must also be set
    % to a value > 0 although it is ignored. Otherwise Auto control is automatically chosen (for backwards compatibility reasons).
    % ziDAQ('set', h, 'bandwidth', 100);
% overwrite_defaults{:,end+1} = {'bandwidth_control'; 2};

    % Sets the bandwidth overlap mode (default 0). If enabled, the bandwidth of a
    % sweep point may overlap with the frequency of neighboring sweep points. The
    % effective bandwidth is only limited by the maximal bandwidth setting and
    % omega suppression. As a result, the bandwidth is independent of the number
    % of sweep points. For frequency response analysis bandwidth overlap should be
    % enabled to achieve maximal sweep speed (default: 0). 0 = Disable, 1 = Enable.
% overwrite_defaults{:,end+1} = {'bandwidth_overlap'; 0};

    %Demodulator time constant.
% overwrite_defaults{:,end+1} = {'demod_time_constant'; 0.007};

    %Demodulator rate.
% overwrite_defaults{:,end+1} = {'demod_rate'; 13e3};

    % Define device channels.
% additional_settings.channels.demod_c = '0'; % demod channel, for paths on the device
% additional_settings.channels.demod_idx = str2double(demod_c)+1; % 1-based indexing, to access the data
% additional_settings.channels.out_c = '0'; % signal output channel
% additional_settings.channels.out_mixer_c = 2; % Define the value of the instrument's default Signal Output mixer channel.
% additional_settings.channels.in_c = '0'; % signal input channel
% additional_settings.channels.osc_c = '0'; % oscillator
% additional_settings.channels.imp_c = '0'; % IA channel
% additional_settings.channels.imp_index = 1; % IA, 1-based indexing, to access the data
end

%% Sweep by
% if strcamp(sweep_param, 'osc_frequency')
% 
% elseif strcamp(sweep_param, 'test_signal')
% 
% elseif strcamp(sweep_param, 'bias_voltage')


% for s = {[sweep_param '_i']}
%     eval(['data.' s '=' s ';'])
% end
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

