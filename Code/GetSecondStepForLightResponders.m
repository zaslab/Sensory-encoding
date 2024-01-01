function [data, data_labels] = GetSecondStepForLightResponders(data, data_labels)
%This function removes the first ON step data of ASH, ASK and ADL 

problem_neurons = ismember(data_labels(:,1), [3, 18, 7, 8, 11, 12]);
first_step = data_labels(:,3) == 1;
second_step = data_labels(:,3) == 3;
data_labels(first_step & problem_neurons, 3) = 3;
data_labels(second_step & problem_neurons, 3) = 1;
data(data_labels(:,3) == 3,:) = [];
data_labels(data_labels(:,3) == 3,:) = [];

end
    