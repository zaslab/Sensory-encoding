function [full_trace, full_ids] = CombineOnOffTraces(data, data_labels)
asf
done = [];
full_trace = [];
full_ids = [];
counter = 1;
for row = 1:(size(data_labels,1) - 1)
    if ismember(row, done)
        continue
    end
    ind = find(sum(data_labels(row,[1,2,4,5,6]) == data_labels(:,[1,2,4,5,6]),2) == 5);

    full_trace(counter,:) = [data(row,:), data(ind(2),:)]; 
    full_ids(counter,:) = data_labels(row,:);

    done = [done, ind(2)];
    counter = counter + 1;
end