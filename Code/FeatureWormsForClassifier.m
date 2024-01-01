% clear
load('relevant_neurons.mat')
load('ordered_neurons.mat')

load('osm6.mat')
load('unc13.mat')
load('unc31.mat')
%%
data = CleanSensoryData(osm6, unc31, unc13);
data = data(1);
[all_neurons_data, variable_matrix] = prepare_PCA_data(data , 'step_benchmark');
%%
var_mat = variable_matrix;
variable_matrix = variable_matrix(:,1:4);
[mean_trace_data, unique_ids] = calculate_step_means(variable_matrix, all_neurons_data, 0,0);
for u = 1:size(all_neurons_data,1)
    all_neurons_data(u,:) = smooth(all_neurons_data(u,:));
end

%%
unique_ids= unique_ids(:,1:4);

new_variable_matrix = variable_matrix;
new_all_neurons_data = all_neurons_data;

data = Reallign(new_all_neurons_data, variable_matrix);
%%
[data, data_labels] = GetSecondStepForLightResponders(data, var_mat);
% [data, data_labels] = CombineOnOffTraces(data, data_labels);
% data_labels(:,3) = 1;

%%
cd('G:\Google Drive\code\catch22-main\wrap_Matlab')
all_trace_features = zeros(size(data,1), 22);
for i = 1:size(data,1)
    all_trace_features(i,:) = catch22_all(smooth(data(i,:)), false);
end

cd('G:\Google Drive\MultiK\All matrices\DataForStatistics')
%%


mean_data = nanmean(all_trace_features, 1);
std_data = nanstd(all_trace_features, [], 1);

% Normalize each time point in each trace
for i = 1:size(all_trace_features,1)
    all_trace_features(i, :) = (all_trace_features(i, :) - mean_data) ./ std_data;
end


%%
 
[coefforth,score,latent,tsquared,explained,mu]= pca(all_trace_features);
[on_diffs, wt_data] = calculate_activations(data_labels, data,0,0, 'all_window', 26:60,'ASH_window', 26:40);




%% Best PCs
[active_variable_matrix,unique_ids1] = CombineNeuronSides(data_labels,unique_ids);
[activations_by_worm, activations_by_worm_identity,all_features_concatenated] = ...
    MakeSyntheticWorms(active_variable_matrix, wt_data, score, 20, [1:3]);
%%
%  save('G:\Google Drive\MultiK\All matrices\DataForStatistics\WithFeatures\osm6synth_activations_by_worm_identity_combos.mat', "activations_by_worm_identity")
%  save('G:\Google Drive\MultiK\All matrices\DataForStatistics\WithFeatures\osm6synth_activations_by_worm_combos.mat', "activations_by_worm")
%  save('G:\Google Drive\MultiK\All matrices\DataForStatistics\WithFeatures\osm6synth_all_PCs_combos.mat', "all_features_concatenated")

%% Worst PCs

[activations_by_worm, activations_by_worm_identity,all_features_concatenated_worst] = ...
    MakeSyntheticWorms(active_variable_matrix, wt_data, score, 20, [20:22]);


% save('osm6synth_all_PCs_combos_worst.mat', "all_features_concatenated_worst")
%%
strain_matrices= {};
bothsteps_matrix = [];
features_to_use = [1:3];
for i = 1:length(features_to_use)

    [all_features,unique_ids]= calculate_combined_trace_features...
        (data_labels, score(:,features_to_use(i)));
    strain_matrices(i,:) = ConstructColorMatrix(unique_ids,all_features, 1);
    % strain_matrices{i,1} = RemoveCrossReads(strain_matrices{i,1}, 1);
    strain_matrices{i,1} = CompleteStrainMatrix(strain_matrices{i,1});
    [matrix, bothstepids ] = MakeWormFeatureVectors(strain_matrices(i,:),0);
    bothsteps_matrix = [bothsteps_matrix, matrix(:,1:13)];

end
%%

bothsteps_matrix = [];
bothstepids = [];
data_labels(:,3) = 1;
features_to_use = [1:3];
for i = 1

    [matrix, bothstepids ] = MakeWormFeatureVectors(strain_matrices(i,:),0);
    bothsteps_matrix = [bothsteps_matrix, matrix(:,1:13)];

end

%%
mean_data = mean(bothsteps_matrix, 1);
std_data = std(bothsteps_matrix, [], 1);

% Normalize each time point in each trace
for i = 1:size(bothsteps_matrix,1)
    bothsteps_matrix(i, :) = (bothsteps_matrix(i, :) - mean_data) ./ std_data;
end
%%
all_PCs= bothsteps_matrix;
% save('G:/Google Drive/MultiK/All matrices/DataForStatistics/all_PCs', 'all_PCs')
% save('G:/Google Drive/MultiK/All matrices/DataForStatistics/bothstepids', 'bothstepids')
%%
function [all_features,unique_ids]= calculate_combined_trace_features...
                                    (variable_matrix, features)

neuron_pairs = [1,nan;
    16,nan;
    2,20;
    3,18;
    4,19;
    5,13;
    6,15;
    7,12;
    8,11;
    9,nan;
    21,nan;
    10,17;
    14,22];
for i = 1:size(neuron_pairs,1)
    if ~isnan(neuron_pairs(i,2))
        variable_matrix(variable_matrix(:,1) == neuron_pairs(i,2)) = neuron_pairs(i,1);
    end
end

counter = 1;

[unique_ids] = unique(variable_matrix(:,1:4), 'rows', 'first');
all_features = {};
for i = 1:size(unique_ids,1)

    neuron_group = unique_ids(i,:);
    this_neuron = neuron_group(1);
    [neuron,~] = find(variable_matrix(:,1) == neuron_group(1));

    [row, col] = find(neuron_pairs == this_neuron);
    if col == 1
        scndneuron = neuron_pairs(row, 2);
    else
        scndneuron = neuron_pairs(row, 1);
    end

    [neuron2,~] = find(variable_matrix(:,1) == scndneuron);


    neuron = cat(1,neuron,neuron2);
    [cond,~] = find(variable_matrix(:,2) == neuron_group(2));
    [step,~] = find(variable_matrix(:,3) == neuron_group(3));
    [strain,~] = find(variable_matrix(:,4) == neuron_group(4));

    point_to_plot = intersect(intersect(intersect(cond,step), strain), neuron);
    wormnums = variable_matrix(point_to_plot, 5);
    sides = variable_matrix(point_to_plot, 6);

    feature = features(point_to_plot);
    

    feature_R = feature(sides == 1);
    feature_L = feature(sides == 2);
    mean_feature = nanmean([feature_L,feature_R],2);
    %


    all_features{counter} = mean_feature;


    counter = counter + 1;
end
all_features = all_features';
end
