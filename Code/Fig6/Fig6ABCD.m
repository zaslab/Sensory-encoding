function Fig6ABCD()

neurons  = {'AWCs','ASJ','ASH','ASG','ADF','ASI','ADL','ASK','ASER','AWA','AWB','AWCw','ASEL'};
cond_names = {'IAA' 'DA' 'Gly' 'NaCl' 'Quin' 'SDS' 'Ctrl'};
set_names = {['_all_' , 'synth_osm6_acts_only', '_all_ind_']
    ['_all_' , 'synth_osm6_PCA_3_features', '_all_ind_']
    ['_all_' , 'synth_osm6_acts_PCA_3_features', '_all_ind_']};

sample_sizes = ones(3,7) * 20;

all_conds_n = [];
all_conds_c = [];

for i = 1:length(set_names)

    set_name = set_names(i);
    which_metric = 2;
    num_repeats = 20;
    all_3_ind_metrics = GetConfusionMetrics(set_name,sample_sizes,'WT',which_metric,num_repeats);

    mean_n_F = mean(all_3_ind_metrics,2);
    
    mean_c_F = permute(mean(all_3_ind_metrics,1), [2,1,3]);
    
    all_conds_n = cat(2, all_conds_n, permute(mean_n_F, [3,2,1]) );
    all_conds_c = cat(2, all_conds_c, permute(mean_c_F, [3,2,1]) );

end

%% Fig 7 A

by_neuron = permute(mean(all_conds_n,1),[2,3,1]);
[~,sorted] = sort(mean(by_neuron ,1),'descend');
sorted_neuron = by_neuron(:,sorted);
figure; imagesc(sorted_neuron')
colormap(othercolor('Blues9', 50))
caxis([min(sorted_neuron(:)), max(sorted_neuron(:))])
caxis([0.3 0.9])
set(gca, 'XTick',[1,2,3] )
set(gca, 'XTickLabel', {'Amplitude', 'Dynamics', 'Both'})
set(gca, 'YTick',[1:13] )
set(gca, 'YTickLabel', neurons(sorted))

%% Fig 7 B
amplitude = by_neuron(1,:)./by_neuron(3,:);
dynamics = by_neuron(2,:)./by_neuron(3,:);

figure;scatter(amplitude,dynamics, 50,'filled')
hold on; plot([0,2], [0,2])

xlim([0.601 1])
ylim([0.6 1])
text(amplitude+0.005,dynamics+0.002,neurons, 'FontSize', 12)
set(gca, 'XTick', 0.5:0.1:1)
set(gca, 'YTick', 0.5:0.1:1)
xlabel('Amplitdue')
ylabel('Dynamics')

%% Fig 7 C
by_stim = permute(mean(all_conds_c,1),[2,3,1]);
[~,sorted] = sort(by_stim(3,:),'descend');
sorted_stim = by_stim(:,sorted);
figure; imagesc(sorted_stim')
set(gca, 'XTick',[1,2,3] )
set(gca, 'XTickLabel', {'Amplitude', 'Dynamics', 'Both'})
set(gca, 'YTickLabel', cond_names(sorted))
colormap(othercolor('Blues9', 50))
caxis([min(sorted_stim(:)), max(sorted_stim(:))])
caxis([0.3 1])

%% Fig 7 D
amplitude = by_stim(1,:)./by_stim(3,:);
dynamics = by_stim(2,:)./by_stim(3,:);

figure;scatter(amplitude,dynamics, 50,'filled')
hold on; plot([0,2], [0,2])

xlim([0.75 1])
ylim([0.75 1])
dynamics(1) = dynamics(1) - 0.01;
amplitude(1) = amplitude(1) - 0.008;
text(amplitude+0.005,dynamics+0.002,cond_names, "FontSize",12)
set(gca, 'XTick', 0.5:0.1:1)
set(gca, 'YTick', 0.5:0.1:1)
xlabel('Amplitdue')
ylabel('Dynamics')
end
%%
function repeat_metrics = GetConfusionMetrics(set_names,sample_size,test_strain,which_metric,num_repeats)

repeat_metrics = [];

for repeat = 1:num_repeats

    all_3_ind_metrics = [];

    for i = 1:size(set_names)
        curr_set = set_names{i};
        confusions = GetConfusions(test_strain, '', curr_set);
        mean_neuron_confusions  = [];

        if iscell(confusions{1})

            for num = 1:length(confusions)
                curr_confusion = confusions{num};
                neuron_confusion = curr_confusion{repeat};
                mean_neuron_confusions = cat(3, mean_neuron_confusions, neuron_confusion);
            end

        else
            mean_neuron_confusions = GetMeanConfusion(confusions);
        end

        all_metrics = GetMetrics(mean_neuron_confusions, sample_size);
        permuted_metrics = permute(all_metrics(which_metric ,:,:), [3,2,1]);
        all_3_ind_metrics = cat(3,all_3_ind_metrics, permuted_metrics);

    end
    repeat_metrics= cat(3, repeat_metrics,all_3_ind_metrics);

end
end

%%
function all_metrics = GetMetrics(mean_neuron_confusions, sample_size)

all_metrics = [];

for j = 0:4:12

    for counter = 1:4

        i= counter + j;

        if i > 13
            continue
        end

        [~, metrics] = ConfusionMetrics(mean_neuron_confusions(:,:,i),sample_size);
        all_metrics = cat(3, all_metrics, metrics);

    end
end


end


%%
function mean_neuron_confusions = GetMeanConfusion(confusions)

mean_neuron_confusions = [];

for num = 1:length(confusions)
    curr_confusions = confusions{num};
    neuron_confusions = [];

    for c = 1:length(curr_confusions)
        neuron_confusions = cat(3, neuron_confusions, curr_confusions);
    end

    mean_neuron_confusions = cat(3, mean_neuron_confusions, mean(neuron_confusions, 3));

end
end
%%
function [weighted_F, metrics] = ConfusionMetrics(cm_test, sample_sizes)

sample_sizes = sample_sizes(:,1:size(cm_test,1));

for i = 1:size(sample_sizes,2)
    cm_test(i,:) = cm_test(i,:) * sample_sizes(i);

end

tp_m = diag(cm_test);

for i = 1:size(sample_sizes,2)

    TP = tp_m(i);
    FP = sum(cm_test(i,:), 2) - TP;
    FN = sum(cm_test(:,i), 1) - TP;
    TN = sum(cm_test(:)) - TP - FP - FN;

    Accuracy(i) = (TP+TN)./(TP+FP+TN+FN);

    TPR(i) = TP./(TP + FN);%tp/actual positive  RECALL SENSITIVITY
    if isnan(TPR(i))
        TPR(i) = 0;
    end
    PPV(i) = TP./ (TP + FP); % tp / predicted positive PRECISION
    if isnan(PPV(i))
        PPV(i) = 0;
    end
    TNR(i) = TN./ (TN+FP); %tn/ actual negative  SPECIFICITY
    if isnan(TNR(i))
        TNR(i) = 0;
    end
    FPR(i) = FP./ (TN+FP);
    if isnan(FPR(i))
        FPR(i) = 0;
    end
    FScore(i) = (2*(PPV(i) .* TPR(i))) ./ (PPV(i)+TPR(i));

    if isnan(FScore(i))
        FScore(i) = 0;
    end
end
metrics = [Accuracy;
    FScore;
    TPR;
    PPV;
    TNR;
    FPR];
weighted_F = sum(sample_sizes .* FScore)/sum(sample_sizes);
end