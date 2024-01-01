clear
load("all_trace_features.mat")
load("traces.mat")
load('all_ids.mat')
load('all_data.mat')

%%
all_trace_features_raw = all_trace_features ;
mean_data = nanmean(all_trace_features, 1);
std_data = nanstd(all_trace_features, [], 1);

% Normalize each time point in each trace
for i = 1:size(all_trace_features,1)
    all_trace_features(i, :) = (all_trace_features(i, :) - mean_data) ./ std_data;
end

[coefforth,score,latent,tsquared,explained,mu]= pca(all_trace_features);

%%
all_trace_features = score;
[unique_ids] = unique(all_ids, 'rows', 'first');
neurons = unique(all_ids(:,1));
num_to_sample = 20;
features_to_use = [1:22];
feat_num = length(features_to_use);
all_features_concatenated = [];
activations_by_worm = [];
activations_by_worm_identity = [];
strain = 1;


for neuron = 1:max(unique_ids(:,1))

    for cond = 1:max(unique_ids(:,2))
        
        rows = ((cond - 1) * num_to_sample + 1): cond*num_to_sample;
        
        for step = 1:2
        
    
    
    [neuron_loc,~] = find(all_ids(:,1) == neuron);
    [cond_loc,~] = find(all_ids(:,2) == cond);
    [step_loc,~] = find(all_ids(:,3) == step);
    
    
    point_to_plot_no_step = intersect(cond_loc,neuron_loc);
    if isempty(point_to_plot_no_step)
        continue
    end
    stim = all_ids(cond_loc(1), 5);
    point_to_plot = intersect(point_to_plot_no_step, step_loc);
    if isempty(point_to_plot)
        continue
    end
    sample_neurons = randsample(length(point_to_plot),num_to_sample, true);
    activities = all_data(point_to_plot(sample_neurons));
    while sum(isnan(activities)> 0)
        nan_num = find(isnan(activities));
        new_samp = randsample(length(point_to_plot),length(nan_num), true);
        activities(nan_num) = all_data(point_to_plot(new_samp));
    end

    activations_by_worm(rows, neuron + (step - 1) * 11) = activities;
    
    activations_by_worm_identity(rows, 1:2) = repmat([stim, cond], num_to_sample,1);
    all_features_concatenated(rows, ((neuron - 1) * feat_num  + 1):(neuron*feat_num ))...
            = all_trace_features(point_to_plot_no_step(sample_neurons),features_to_use);
        end
    end
end

%%
all_features_concatenated(sum(activations_by_worm == 0, 2) > 10,:) = [];
activations_by_worm_identity(sum(activations_by_worm == 0, 2) > 10,:) = [];
activations_by_worm(sum(activations_by_worm == 0, 2) > 10,:) = [];

%%
row_mean = mean(traces(:,21:end),2);
for i = 1:size(traces(:,21:end),2)
    traces(:, i) = (traces(:, i) - row_mean );
end

mean_data = mean(traces, 1);
std_data = std(traces, [], 1);

% Normalize each time point in each trace
for i = 1:size(traces,1)
    traces(i, :) = (traces(i, :) - mean_data) ./ std_data;
end

