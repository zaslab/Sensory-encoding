function [sortedData] = sortAWCs(data)

for i = 1:size(data,3)
    if sum(isnan(data(1,:,i))) > 1 || sum(isnan(data(16,:,i))) > 1
        data(1,:,i) = NaN;
        data(16,:,i) = NaN;
    else if abs(nanmean(data(1,:,i))) < abs(nanmean(data(16,:,i)))
            data(:,:,i) = data([16,2:15,1,17:22],:,i);
        end
    end
end
sortedData = data;
end