function [sweep_range, sweep_pts, frequency_vec, amplitude_vec, offset_vec]  = MFIA_freq_amp_bias_value_pairs_withParser(sweep_order, varargin)


frequency_vec = [];
amplitude_vec = [];
offset_vec = [];

% Define parameters relevant to this example. Default values specified by the
% inputParser below are overwritten if specified as name-value pairs via the
% `varargin` input argument.
p = inputParser;
p.KeepUnmatched=true;
isnonnegscalar = @(x) isnumeric(x) && isscalar(x) && (x > 0);

p.addParameter('start_frequency', 1e2, @isnumeric);
p.addParameter('stop_frequency', 500e3, @isnumeric);
p.addParameter('pts_frequency', 100, isnonnegscalar);
p.addParameter('freq_xmapping', 1, isnonnegscalar);

p.addParameter('start_amplitude', 0.05, @isnumeric);
p.addParameter('stop_amplitude', 0.25, @isnumeric);
p.addParameter('pts_amplitude', 21, isnonnegscalar);

p.addParameter('start_offset', 0, @isnumeric);
p.addParameter('stop_offset', -1, @isnumeric);
p.addParameter('pts_offset', 21, isnonnegscalar);


p.parse(varargin{:});
sw1 = sweep_order{1};
eval(['sweep_range=[p.Results.start_' sw1 ' p.Results.stop_' sw1 '];'])
eval(['sweep_pts=p.Results.pts_' sw1 ';'])

start1 = eval(['p.Results.start_' sw1]);
stop1 = eval(['p.Results.stop_' sw1]);
pts1 = eval(['p.Results.pts_' sw1]);
vec1 = linspace(start1,stop1,pts1);
if strcmpi(sw1, 'frequency') && p.Results.freq_xmapping
    vec1 = logspace(log10(start1),log10(stop1),pts1);
end


sw2 = sweep_order{2};
start2 = eval(['p.Results.start_' sw2]);
stop2 = eval(['p.Results.stop_' sw2]);
pts2 = eval(['p.Results.pts_' sw2]);
vec2 = linspace(start2,stop2,pts2);
if strcmpi(sw2, 'frequency') && p.Results.freq_xmapping
    vec2 = logspace(log10(start2),log10(stop2),pts2);
end


sw3 = sweep_order{3};
start3 = eval(['p.Results.start_' sw3]);
stop3 = eval(['p.Results.stop_' sw3]);
pts3 = eval(['p.Results.pts_' sw3]);   
vec3 = linspace(start3,stop3,pts3);
if strcmpi(sw3, 'frequency') && p.Results.freq_xmapping
    vec3 = logspace(log10(start3),log10(stop3),pts3);
end

[mat3, mat2] = meshgrid(vec3,vec2);
eval([sw1 '_vec=vec1'';'])
eval([sw2 '_vec=mat2(:);'])
eval([sw3 '_vec=mat3(:);'])
end

