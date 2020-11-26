function [sweep_range, sweep_pts, frequency_vec, amplitude_vec, offset_vec]  = MFIA_freq_amp_bias_value_pairs_withparser(sweep_order, varargin)

% Define parameters relevant to this example. Default values specified by the
% inputParser below are overwritten if specified as name-value pairs via the
% `varargin` input argument.
p = inputParser;
p.KeepUnmatched=true;
isnonnegscalar = @(x) isnumeric(x) && isscalar(x) && (x > 0);

p.addParameter('start_frequency', 1e2, @isnumeric);

stop = 500e3;
% if strfind(props.devicetype, 'MF')
%     stop = 10e6;
% end
p.addParameter('stop_frequency', stop, @isnumeric);
p.addParameter('pts_frequency', 100, isnonnegscalar);

p.addParameter('start_amplitude', 0.05, @isnumeric);
p.addParameter('stop_amplitude', 0.3, @isnumeric);
p.addParameter('pts_amplitude', 20, isnonnegscalar);

p.addParameter('start_offset', 1, @isnumeric);
p.addParameter('stop_offset', -1, @isnumeric);
p.addParameter('pts_offset', 20, isnonnegscalar);


p.parse(varargin{:});
sw1 = sweep_order{1};
eval(['sweep_range=[p.Results.start_' sw1 ' p.Results.stop_' sw1 '];'])
eval(['sweep_pts=p.Results.pts_' sw1 ';'])


sw2 = sweep_order{2};
start2 = eval(['start_' sw2]);
stop2 = eval(['stop_' sw2]);
pts2 = eval(['pts_' sw2]);
vec2 = linspace(start2,stop2,pts2);  

sw3 = sweep_order{3};
start3 = eval(['start_' sw3]);
stop3 = eval(['stop_' sw3]);
pts3 = eval(['pts_' sw3]);   
vec3 = linspace(start3,stop3,pts3);

[mat3, mat2] = meshgrid(vec2,vec3);
eval([sw2 '_vec=' mat2(:) ';'])
eval([sw3 '_vec=' mat3(:) ';'])
end

