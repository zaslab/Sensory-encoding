function [all_act_diffs, unordered] = calculate_activations(...
    variable_matrix, all_neurons_data, combine_neurons, varargin)
% Calculate activation magnitude for each neuron
% Neuron identity recoded in unique_ids

if isempty(varargin) 
    window = 26:58;
else
    for i = 1:2:length(varargin)
        if strcmp('all_window', varargin{i})
            all_window = varargin{i+1};
        end
        if strcmp('ASH_window', varargin{i})
            ASH_window = varargin{i+1};
        end
    end
end

unordered = zeros(size(variable_matrix,1),1);
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

[unique_ids] = unique(variable_matrix, 'rows', 'first');
all_act_diffs = {};
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
    
    isASHoff = (neuron_group(1)== 18 | neuron_group(1)== 3) & neuron_group(3)== 2;
    if isASHoff
        window = ASH_window;
    else 
        window = all_window;
    end
    pre_window = (all_window(1) - all_window(1) + 1):(all_window(1) - 5);
    [~, max_ind] = max(all_neurons_data(point_to_plot,window ),[],2);
    post_act_mean = mean(all_neurons_data(point_to_plot,window ),2);
    [~, min_ind] = min(all_neurons_data(point_to_plot,window ),[],2);
    pre_act = mean(all_neurons_data(point_to_plot,pre_window),2);
    min_ind = min_ind + window(1);
    max_ind = max_ind + window(1);
    post_act = [];
    

    for worm = 1:length(post_act_mean)
        if nanmean(post_act_mean) < nanmean(pre_act) & ~isASHoff
            post_act(worm) = mean(all_neurons_data(point_to_plot(worm),(min_ind(worm) - 1):(min_ind(worm) + 1)),2);
           
        else
            post_act(worm) = mean(all_neurons_data(point_to_plot(worm),(max_ind(worm) - 1):(max_ind(worm) + 1)),2);
           
        end
    end
    
    act_diff = post_act' - pre_act;

    unordered(point_to_plot) = act_diff;
    
 
    all_act_diffs{counter} = act_diff;
    counter = counter + 1;
end
all_act_diffs = all_act_diffs';
end
