function cond_loc = FindCondNum(data, current_cond)

where_cond = cellfun(@(x) contains(x, current_cond), data(1:7,3), 'UniformOutput', false);
cond_locs = find(cell2mat(where_cond));
if length(cond_locs) > 1 
    for i = 1:length(cond_locs)
        if strcmp(current_cond, data(cond_locs(i),3))
            cond_loc = cond_locs(i);
        end
    end
else
    cond_loc = cond_locs;
end
end