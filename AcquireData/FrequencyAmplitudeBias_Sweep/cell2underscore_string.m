function out = cell2underscore_string(sw)
out = [];
for c = sw
    out = [out c{:} '_'];
end
out(end) = [];
end