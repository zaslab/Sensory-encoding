function [activations_by_worm, activations_by_worm_identity,all_features_concatenated]...
          = MakeSyntheticWorms(all_ids, all_data, all_trace_features,...
            num_to_sample,features_to_use)

[unique_ids] = unique(all_ids, 'rows', 'first');

asf
feat_num = length(features_to_use);
on_features_concatenated = [];
off_features_concatenated = [];
activations_by_worm = [];
activations_by_worm_identity = [];
relevant_neurons = unique(unique_ids(:,1));
% plotcounter = 1;
% figure
for neuron = 1:length(relevant_neurons)
    neuron_id = relevant_neurons(neuron);
    feat_cols = ((neuron - 1) * feat_num  + 1):(neuron*feat_num );
    for cond = 1:max(unique_ids(:,2))
        
        
        rows = ((cond - 1) * num_to_sample + 1): cond*num_to_sample;
        
        for step = 1:2
        
    
    
    [neuron_loc,~] = find(all_ids(:,1) == neuron_id );
    [cond_loc,~] = find(all_ids(:,2) == cond);
    [step_loc,~] = find(all_ids(:,3) == step);
    
    
    point_to_plot_no_step = intersect(cond_loc,neuron_loc);
    if isempty(point_to_plot_no_step)
        continue
    end

    point_to_plot = intersect(point_to_plot_no_step, step_loc);
    if isempty(point_to_plot)
        continue
    end
    sample_neurons = randsample(length(point_to_plot),num_to_sample, true);
    activities = all_data(point_to_plot(sample_neurons));
    if step == 1
    on_features_concatenated(rows,feat_cols)...
            = all_trace_features(point_to_plot(sample_neurons),features_to_use);
    else
    off_features_concatenated(rows,feat_cols)...
            = all_trace_features(point_to_plot(sample_neurons),features_to_use);

    end

    while sum(isnan(activities)> 0)
        nan_num = find(isnan(activities));
        new_samp = randsample(length(point_to_plot),length(nan_num), true);
        activities(nan_num) = all_data(point_to_plot(new_samp));
        if step == 1
        on_features_concatenated(nan_num + rows(1) - 1,feat_cols) = ...
            all_trace_features(point_to_plot(new_samp),features_to_use);
        
        else
            off_features_concatenated(nan_num + rows(1) - 1,feat_cols) = ...
            all_trace_features(point_to_plot(new_samp),features_to_use);

    end
    end

    activations_by_worm(rows, neuron + (step - 1) * 13) = activities;
    
    activations_by_worm_identity(rows, 1) = repmat( cond, num_to_sample,1);
%     half_len = length(point_to_plot_no_step)/2;
%     point_to_plot_no_step((half_len+1):end) = point_to_plot_no_step(1:half_len);
    
    %     if plotcounter < 101
%     subplot(10,10,plotcounter)
%     histogram(all_data(point_to_plot),'NumBins',5)
%     plotcounter = plotcounter + 1;
%     end
        end
    end
end
all_features_concatenated = cat(2, on_features_concatenated, off_features_concatenated);
end