function out = compress_array_exp(array,a,b,varargin)
Tflag = false;
if size(array,2)>size(array,1)
    array = array.';
    Tflag = true;
end
p = inputParser;
p.addParameter('Reverse', false);
p.addParameter('Interp', false);
p.addParameter('Method', 'linear', @isstring);
p.parse(varargin{:});
ind = unique(round((a.^linspace(0,b,length(array))-1)*length(array)/(a^b-1)));
if p.Results.Reverse
    ind = fliplr(abs(ind-length(array)));
end
out = array(ind(ind>0),:);
if p.Results.Interp
    out = interp1(out, 1:(length(out)/length(array)):length(out), p.Results.Method);
end
if Tflag || size(out,2)==1
    out = out.';
end