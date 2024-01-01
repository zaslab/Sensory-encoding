function [mean_trace_data, unique_ids] = calculate_step_means(variable_matrix, all_neurons_data, combine_neurons)
% Calculate mean trace for each neuron in each cond.
% Neuron identity recorded in unique_ids

% combine_neurons = 1 for combined L and R. Also combines ASEs and AWCs
% which is bad. Currently spits out double results (twice for each pair)

neuron_pairs = [1,16;
    2,20;
    3,18;
    4,19;
    5,13;
    6,15;
    7,12;
    8,11;
    9,21;
    10,17;
    14,22];

counter = 1;
mean_trace = [];
[unique_ids] = unique(variable_matrix, 'rows', 'first');

for i = 1:size(unique_ids,1)

    neuron_group = unique_ids(i,:);
    this_neuron = neuron_group(1);
    [neuron,~] = find(variable_matrix(:,1) == neuron_group(1));
    if combine_neurons
        [row, col] = find(neuron_pairs == this_neuron);
        if col == 1
            scndneuron = neuron_pairs(row, 2);
        else
            scndneuron = neuron_pairs(row, 1);
        end

        [neuron2,~] = find(variable_matrix(:,1) == scndneuron);


        neuron = cat(1,neuron,neuron2);

    end
    [cond,~] = find(variable_matrix(:,2) == neuron_group(2));

    [step,~] = find(variable_matrix(:,3) == neuron_group(3));
    [strain,~] = find(variable_matrix(:,4) == neuron_group(4));
    point_to_plot = intersect(intersect(intersect(cond,step), strain), neuron);


    PCgroup = all_neurons_data(point_to_plot ,:);
    for i = 1:size(PCgroup,1)
        PCgroup(i,:) = PCgroup(i,:) - 1;
    end

    current_mean_trace = nanmean(PCgroup,1);
    current_std_trace = nanstd(PCgroup,1);
    group_size = sum(~isnan(PCgroup(:,1)));
    current_ste_trace = current_std_trace./sqrt(group_size);

    mean_trace(counter,:) = current_mean_trace;
    std_trace(counter,:) = current_std_trace;
    ste_trace(counter,:) = current_ste_trace;


    counter = counter + 1;
end
mean_trace_data = {mean_trace, std_trace, ste_trace};
end