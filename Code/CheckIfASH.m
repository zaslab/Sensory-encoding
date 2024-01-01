function aligned_trace = CheckIfASH(data, max_change)

[peak_size, peak] = max(data(30:50));
start = mean(data(1:10));
post_step = mean(data((max_change + 30):(max_change + 40)));
if  peak_size > (start + 0.08) & post_step < (start - 0.1)
    max_change = max([peak - 8,1]);

end
aligned_trace = data((max_change) : (max_change + 139));
end
