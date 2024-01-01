asf
neurons = unique(all_ids(:,1));

num_to_sample = 50;
% features_to_use = [4,8,10,12,3,17];
features_to_use = [4,8,10];
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
    [strain_loc,~] = find(all_ids(:,4) == strain);
    
    point_to_plot_no_step = intersect(intersect(cond_loc,neuron_loc), strain_loc);

    point_to_plot = intersect(point_to_plot_no_step, step_loc);
    sample_neurons = randsample(length(point_to_plot),num_to_sample, true);
    activities = all_data(point_to_plot(sample_neurons));
    while sum(isnan(activities)> 0)
        nan_num = find(isnan(activities));
        new_samp = randsample(length(point_to_plot),length(nan_num), true);
        activities(nan_num) = all_data(point_to_plot(new_samp));
    end
    activations_by_worm(rows, neuron + (step - 1) * 11) = activities;
    
    activations_by_worm_identity(rows, 1:2) = repmat([neuron, cond], num_to_sample,1);
    all_features_concatenated(rows, ((neuron - 1) * feat_num  + 1):(neuron*feat_num ))...
            = all_trace_features(point_to_plot_no_step(sample_neurons),features_to_use);
        end
    end
end

%%


save( 'aravi_all_features_concatenated','all_features_concatenated')
save( 'aravi_activations_by_worm','activations_by_worm')
save('aravi_activations_by_worm_identity', 'activations_by_worm_identity')
