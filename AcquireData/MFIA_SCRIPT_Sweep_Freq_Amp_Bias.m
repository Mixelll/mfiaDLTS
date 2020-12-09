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

% Sample name
sample_name = 'B5 b5 150um 14 C';

% Save path
save_path = 'C:\Users\Public\Documents\MATLAB\mfiaDLTS\AcquireData\vars';
desired_order = {'offset','frequency', 'amplitude'};

% create slice planes 
slice_planes = {};
slice_planes(:,end+1) = {'frequency'; [100 10000 500e3]};
slice_planes(:,end+1) = {'amplitude'; []}; 
slice_planes(:,end+1) = {'offset'; [-0.5]}; 

% Selected data to read (grid = sweep parameter)
% read_param_struct.demod.grid = true;
read_param_struct.demod.r = false;
read_param_struct.demod.phase = false;
read_param_struct.impedance.grid = true;
read_param_struct.impedance.param0 = true;
read_param_struct.impedance.param1 = true;

% Sweep by Frequency and iterate over AC amplitude and offset
sweep_order = {'frequency','amplitude','offset'};
freq_xmapping = {'freq_xmapping', 1}; % set 0 for linear distribution between start and stop, set 1 for log distribution
start_frequency = 100; stop_frequency = 500e3; pts_frequency = 100; % Hz
start_amplitude = 0.05; stop_amplitude = 0.2; pts_amplitude = 20; % V
start_offset = 0; stop_offset = -1; pts_offset = 20; % V

plt_log_freq = true; % set to plot sweep/desired order with logarithmic frequency values
plot_desired_order = true;
plot_sweep_order = true;
plt_cmds = {};
% plt_cmds{end+1} = 'grid on';
plt_cmds{end+1} = 'colorbar(''eastoutside'')';
%% Overwite Defaults (uncomment and change value)

overwrite_defaults = {}; % don't touch
additional_settings = struct; % don't touch

    % BY SETTING THIS, THE INTERNAL SCRIPT DEFAULTS WILL OVERWRITE LABONE GUI INPUT
    % THE INTERNAL DEFAULTS ARE INSIDE MFIA_freq_amp_bias_sweep.
    % VALUES YOU UN-COMMENT BELOW WILL OVERWRITE REGARDLESS of "enable_default"
additional_settings.enable_default = true;
%% Graph and text display settings
    % Graphs
% additional_settings.display.graph.disp = true;
% additional_settings.display.graph.during_sweep = true;
    % Text
% additional_settings.display.text.major.disp = true;
% additional_settings.display.text.major.each_sweep = true;
% additional_settings.display.text.minor.disp = true;
% additional_settings.display.text.minor.each_sweep = false;
%% Saving the plot of the sweeper output (can be hundreds of images)
additional_settings.display.graph.save.if = true;
mkdir(['C:\Users\Public\Documents\MATLAB\mfiaDLTS\AcquireData\vars\' sample_name])
additional_settings.display.graph.save.path = ['C:\Users\Public\Documents\MATLAB\mfiaDLTS\AcquireData\vars\' sample_name];
%% 
%% MF and IA settings
    % IA precision -> measurement speed: 0 - low->fast, 1 - high->medium,
    % 2 - very high->slow
% overwrite_defaults(:,end+1) = {'IA_precision'; 1};

    % Set IA parameter extraction model (from impedance). 
    % 0 - Rp Cp, 1 - Rs Cs, 2 - Rs Ls, 3 - G B, 4 - D Cs,
    % 5 -  Qs Cs, 6 - D Ls, 7 - Q Ls, 8 - Rp Lp, 9 - D Cp
% overwrite_defaults(:,end+1) = {'model'; 0};

    % Enable two terminal measurement.
% overwrite_defaults(:,end+1) = {'two_terminal'; 1};
    
    % Enable two terminal cable length compensation.
% overwrite_defaults(:,end+1) = {'cable_length'; 1};

    % Enable high-pass filter.
% overwrite_defaults(:,end+1) = {'AC'; 0};

    % Enable 50ohm output impedance. Disabled state is 10M ohm.
% overwrite_defaults(:,end+1) = {'imp50ohm'; 0};

    % Enable auto range.
% overwrite_defaults(:,end+1) = {'auto_range'; 0};

    % Current range.
overwrite_defaults(:,end+1) = {'current_range'; 1e-3};

    % Voltage range.
% overwrite_defaults(:,end+1) = {'voltage_range'; 3};

    % Demodulator time constant.
% overwrite_defaults(:,end+1) = {'demod_time_constant'; 0.007};

    % demod data transfer rate, [Hz].
% overwrite_defaults(:,end+1) = {'demod_rate'; 13.39e3};

    % For IA:  Selects the filter roll off to use for the sweep in fixed bandwidth mode. Range between 6 dB/oct and 48 dB/ oct.
% overwrite_defaults(:,end+1) = {'demod_LFP_order'; 8};
%% Define device channels.
% additional_settings.channels.demod_c = '0'; % demod channel, for paths on the device
% additional_settings.channels.demod_idx = str2double(demod_c)+1; % 1-based indexing, to access the data
% additional_settings.channels.out_c = '0'; % signal output channel
% additional_settings.channels.out_mixer_c = 2; % Define the value of the instrument's default Signal Output mixer channel.
% additional_settings.channels.in_c = '0'; % signal input channel
% additional_settings.channels.osc_c = '0'; % oscillator
% additional_settings.channels.imp_c = '0'; % IA channel
% additional_settings.channels.imp_index = 1; % IA, 1-based indexing, to access the data
%% Sweeper settings
    % Sweep timeout.
% overwrite_defaults(:,end+1) = {'timeout'; 120};


    % Perform one single sweep.
% overwrite_defaults(:,end+1) = {'loopcount'; 1};

    % Logarithmic sweep mode.
% overwrite_defaults(:,end+1) = {'xmapping'; 0};

    % Binary scan type.
% overwrite_defaults(:,end+1) = {'scan'; 1};

    % Minimum wait time in seconds between a sweep parameter change and the recording of the next sweep point. This
    % parameter can be used to define the required settling time of the experimental setup. The effective wait time
    % is the maximum of this value and the demodulator filter settling time determined from the Inaccuracy value specified
% overwrite_defaults(:,end+1) = {'settling_time'; 0};

    % Match sweeper accuracy and averaging settings to same settings in the IA
sweep_precision = struct('sweep_inaccuracy',0.01, 'averaging_sample',20, 'averaging_time',0.1, ... % don't touch
'averaging_time_constant',15, 'bandwidth',10, 'max_bandwidth',100, 'omega_suppression', 80); % don't touch
if any(strcmpi(overwrite_defaults, 'IA_precision'))
    switch overwrite_defaults{2, strcmpi(overwrite_defaults, 'IA_precision')}
    case 0
        sweep_precision.sweep_inaccuracy = 0.01;
        sweep_precision.averaging_sample = 20;
        sweep_precision.averaging_time = 0.01;
        sweep_precision.averaging_time_constant = 5;
        sweep_precision.bandwidth = 100;
        sweep_precision.max_bandwidth = 1000;
        sweep_precision.omega_suppression = 60;
    case 1
        sweep_precision.sweep_inaccuracy = 0.01;
        sweep_precision.averaging_sample = 20;
        sweep_precision.averaging_time = 0.1;
        sweep_precision.averaging_time_constant = 15;
        sweep_precision.bandwidth = 10;
        sweep_precision.max_bandwidth = 100;
        sweep_precision.omega_suppression = 80;
    case 2
        sweep_precision.sweep_inaccuracy = 0.0001;
        sweep_precision.averaging_sample = 20;
        sweep_precision.averaging_time = 1;
        sweep_precision.averaging_time_constant = 25;
        sweep_precision.bandwidth = 1;
        sweep_precision.max_bandwidth = 10;
        sweep_precision.omega_suppression = 120;
    end
end

    % overwrite sweep precision
% sweep_precision.sweep_inaccuracy = 0.01;
% sweep_precision.averaging_sample = 20;
% sweep_precision.averaging_time = 0.1;
% sweep_precision.averaging_time_constant = 15;
% sweep_precision.bandwidth = 10;
% sweep_precision.max_bandwidth = 100;
% sweep_precision.omega_suppression = 80;

    % Demodulator filter settling inaccuracy defining the wait time between a sweep parameter change and
    % recording of the next sweep point. The settling time is calculated as the time required to attain the specified
    % remaining proportion [1e-13,0.1] of an incoming step function. Typical inaccuracy the number of filter time
    % constants the sweeper has to wait. The maximum between this value and the settling time is taken as wait time until the
    % next sweep point is recorded values: 10 m for highest sweep speed for large signals, 100 u for precise amplitude
    % measurements, 100 n for precise noise measurements. Depending on the order the settling accuracy will define
%overwrite_defaults(:,end+1) = {'sweep_inaccuracy'; sweep_precision.sweep_inaccuracy};

    % Sets the number of data samples per sweeper parameter point that is considered in the measurement. The maximum
    % between samples, time and number of time constants is taken as effective calculation time.
% overwrite_defaults(:,end+1) = {'averaging_samples'; sweep_precision.averaging_sample};

    % Sets the time during which data samples are processed. The maximum between samples, time and number
    % of time constants is taken as effective calculation time
% overwrite_defaults(:,end+1) = {'averaging_time'; sweep_precision.averaging_time};

    % Sets the effective measurement time per sweeper parameter point that is considered in the
    % measurement. The maximum between samples, time and number of time constants is
    % taken as effective calculation time.
% overwrite_defaults(:,end+1) = {'averaging_time_constant'; sweep_precision.averaging_time_constant};

    % Automatically is recommended in particular for logarithmic sweeps and assures the whole spectrum is covered.
    % Auto: All bandwidth settings of the chosen demodulators are automatically adjusted. For logarithmic sweeps the
    % measurement bandwidth is adjusted throughout the measurement.
    % Fixed: Define a certain bandwidth which is taken for all chosen demodulators for the course of the measurement.
    % Manual: The sweeper module leaves the demodulator bandwidth settings entirely untouched.
% overwrite_defaults(:,end+1) = {'bandwidth_control'; 2};

    % If enabled the bandwidth of a sweep point may overlap with the frequency of neighboring sweep points.
    % The effective bandwidth is only limited by the maximal bandwidth setting and omega suppression. As a result, the bandwidth is independent of
    % the number of sweep points. For frequency response analysis bandwidth overlap should be enabled to achieve maximal sweep speed.
% overwrite_defaults(:,end+1) = {'bandwidth_overlap'; 0};

    % NEP [Hz] Defines the measurement bandwidth for Fixed bandwidth sweep mode, and corresponds to either noise
    % equivalent power bandwidth (NEP), time constant (TC) or 3 dB bandwidth (3 dB) depending on selection. 
% overwrite_defaults(:,end+1) = {'bandwidth'; sweep_precision.bandwidth};

    % [Hz] Limit of the maximum bandwidth used on the demodulator filter. Values above 1 kHz can heavily
    % diminish measurement accuracy in the highfrequency region where the amplitude is no more constant over frequency.
% overwrite_defaults(:,end+1) = {'max_bandwidth'; sweep_precision.max_bandwidth};

    % [dB] Suppression of the omega and 2-omega components. Small omega suppression can diminish measurements of
    % very low or high impedance because the DC component can become dominant. Large omega suppression will have
    % a significant impact on sweep time especially for low filter orders.
% overwrite_defaults(:,end+1) = {'omega_suppression'; sweep_precision.omega_suppression};

    % For sweeper: Selects the filter roll off to use for the sweep in fixed bandwidth mode. Range between 6 dB/oct and 48 dB/ oct.
% overwrite_defaults(:,end+1) = {'sweep_LFP_order'; 8};

%% Run Measurement extract and reshape data
% override defaults set in MFIA_freq_amp_bias_value_pairs_withParser. 
[sweep_range, sweep_pts, frequency_vec, amplitude_vec, offset_vec] = MFIA_freq_amp_bias_value_pairs_withParser(sweep_order, 'start_frequency', start_frequency,...
    'stop_frequency', stop_frequency, 'pts_frequency', pts_frequency, 'start_amplitude', start_amplitude, 'stop_amplitude', stop_amplitude,...
    'pts_amplitude', pts_amplitude, 'start_offset', start_offset, 'stop_offset', stop_offset, 'pts_offset', pts_offset, freq_xmapping{:});

[select_data_sweep_order_struct_vec, full_data_sweep_order_struct_vec] = MFIA_freq_amp_bias_sweep(device_id, additional_settings, sweep_order, sweep_range, sweep_pts, frequency_vec, amplitude_vec, offset_vec, read_param_struct, overwrite_defaults{:});

[select_data_desired_order_3D, select_data_sweep_order_3D] = MFIA_data_reshape_3D(select_data_sweep_order_struct_vec, desired_order, sweep_order, pts_frequency, pts_amplitude, pts_offset, frequency_vec, amplitude_vec, offset_vec);

save([save_path '\' sample_name '_struct_3D_desired_order_' sweep_order_string(desired_order)], 'select_data_desired_order_3D');
save([save_path '\' sample_name '_struct_3D_sweep_order_' sweep_order_string(sweep_order)], 'select_data_sweep_order_3D');
save([save_path '\' sample_name '_struct_vec_sweep_order_' sweep_order_string(sweep_order)], 'select_data_sweep_order_struct_vec');

if plt_log_freq
    desired_order{contains(desired_order, 'frequency')} = 'log_frequency';
    sweep_order{contains(sweep_order, 'frequency')} = 'log_frequency';
    slice_planes(:,contains(slice_planes(1,:), 'frequency')) = {'log_frequency'; log10(slice_planes{2,contains(slice_planes(1,:), 'frequency')})};
end
    
if plot_desired_order
    plot_data3D(select_data_desired_order_3D, desired_order, slice_planes, plt_cmds);
end
if plot_sweep_order
    plot_data3D(select_data_sweep_order_3D, sweep_order, slice_planes, plt_cmds);
end
%%
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

function out = sweep_order_string(sw)
out = [];
for c = sw
    out = [out c{:} '_'];
end
out(end) = [];
end

% Plot data
function plot_data3D(struct, order, slic, plt_cmds)
data = struct.data;
axes = struct.axes;
data_cell = fn_struct2cell(data);
axes_cell = fn_struct2cell(axes);
% Frequency, amplitude and offset grid 
X = axes_cell{2,strcmpi(axes_cell(4,:), order{1})};
Y = axes_cell{2,strcmpi(axes_cell(4,:), order{2})};
Z = axes_cell{2,strcmpi(axes_cell(4,:), order{3})};
% Axis labels
units = {'frequency','log_frequency' 'amplitude', 'offset';' [Hz]', ' [log(Hz)]', ' [V]', ' [V]'};
Xlbl = [strrep(order{1},'_','\_') units{2,strcmpi(units(1,:), order{1})}];
Ylbl = [strrep(order{2},'_','\_') units{2,strcmpi(units(1,:), order{2})}];
Zlbl = [strrep(order{3},'_','\_') units{2,strcmpi(units(1,:), order{3})}];
% Slice planes
Xsl = slic{2,strcmpi(slic(1,:), order{1})};
Ysl = slic{2,strcmpi(slic(1,:), order{2})};
Zsl = slic{2,strcmpi(slic(1,:), order{3})};

figure
clf
[subrow, subcol] = subplot_min_rectangle(size(data_cell,2));
for i = 1:size(data_cell,2)
    subplot(subrow, subcol, i)
    slice(Y,X,Z, data_cell{2,i}, Ysl, Xsl, Zsl);
    xlabel(Ylbl)
    ylabel(Xlbl)
    zlabel(Zlbl)
    title(data_cell{4,i})
    for c = plt_cmds
        eval(c{:});
    end
end
end
