function aligned = Reallign(data, variables, varargin)
% This function aligns the traces to peak activation for more accurate
% calculation of change in activity caused by the step.

normalize = 0;
for i=1:1:length(varargin)
    if strcmp(varargin{i},'normalize') 
        normalize = 1;
    end
end
trace_diffs = diff(data(:,30:50), 1,2); 
aligned = zeros(size(data,1), 140);
[~,max_change] = max(abs(trace_diffs),[],2);
keep_track = [];
for i = 1:length(max_change)


if (variables(i,1) == 3 || variables(i,1) == 18) & variables(i, 3) == 2
    aligned_trace = CheckIfASH(data(i,:), max_change(i));
    keep_track(end + 1) = i;

else 
    aligned_trace = data(i,(max_change(i)) : (max_change(i) + 139));
end

if normalize
aligned_trace = aligned_trace - mean(aligned_trace(20:24));
end

aligned(i,:) = aligned_trace; 

end
end
