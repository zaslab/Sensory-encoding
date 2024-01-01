%%
function [mean_trace,mean_trace_in_PC,unique_ids] = CalculateMeanActiveTrace(data,all_active_neurons,active_variable_matrix,coefforth,combine)
% This function calculates the mean trace of each neuron in each stimulus
% and projects the trace onto PC space
% Use whoacts for only active neurons
counter = 1;
load("relevant_neurons.mat")
load('whoacts.mat')

for i = 1:size(whoacts,1)

    curr_cond = whoacts{i,1};
    curr_neur = whoacts{i,2};
    curr_step = whoacts{i,3};
    cond_num = FindCondNum(data, curr_cond);
    neuron_num = FindNeuronNum(curr_neur, relevant_neurons);
    for k = 1:length(neuron_num)
        toinclude(counter,:) = {neuron_num(k), cond_num, curr_step, 1};
        if max(active_variable_matrix(:,4)) > 1
            toinclude(counter + 1,:) = {neuron_num(k), cond_num, curr_step, 2};
            toinclude(counter + 2,:) = {neuron_num(k), cond_num, curr_step, 3};
            counter = counter + 3;
        else
            counter = counter + 1;
        end
    end
end

if combine == 1
    [active_variable_matrix,toinclude_mat] = CombineNeuronSides(active_variable_matrix,cell2mat(toinclude));
    toinclude_cell = {};
    for i = 1:size(toinclude_mat,1)
        toinclude_cell{i,1} = toinclude_mat(i,1);
        toinclude_cell{i,2} = toinclude_mat(i,2);
        toinclude_cell{i,3} = toinclude_mat(i,3);
        toinclude_cell{i,4} = toinclude_mat(i,4);
    end
    toinclude = toinclude_cell;
end

counter = 1;
mean_trace_in_PC = [];
mean_trace = [];
for i = 1:size(toinclude,1)
    neuron_group = toinclude(i,:);
    [neuron,~] = find(active_variable_matrix(:,1) == neuron_group{1});
    [cond,~] = find(active_variable_matrix(:,2) == neuron_group{2});


    [step,~] = find(active_variable_matrix(:,3) == neuron_group{3});
    [strain,~] = find(active_variable_matrix(:,4) == neuron_group{4});
    point_to_plot = intersect(intersect(intersect(cond,step), strain), neuron);


    PCgroup = all_active_neurons(point_to_plot ,:);

    current_mean_trace = mean(PCgroup,1);
    if size(current_mean_trace,2) == size(coefforth,2)
        mean_trace_in_PC(counter,:) = current_mean_trace * coefforth;
    else
        mean_trace_in_PC=[];
    end
    mean_trace(counter,:) = current_mean_trace;
    counter = counter + 1;
end

unique_ids = cell2mat(toinclude);
end
