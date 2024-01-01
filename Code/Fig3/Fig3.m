%%
function Fig3(data)


for i = 1:7
    if i == 4
        step = 1;
    else
        step = 2;
    end
    data(i,:) = sortByStrength(data(i,:), 1,{'AWC'},step);
end

[all_neurons_data, variable_matrix] = prepare_PCA_data(data, 'step_benchmark',1);
% [mean_trace_data, unique_ids] = calculate_step_means(variable_matrix, all_neurons_data, 0);
for u = 1:size(all_neurons_data,1)
    all_neurons_data(u,:) = smooth(all_neurons_data(u,:));
end

%%
% Use 2nd ON step for light-responding neurons
active_variable_matrix = variable_matrix;
all_neurons_data(active_variable_matrix(:,3) == 1 & ismember(active_variable_matrix(:,1), [3, 18, 7, 8, 11, 12]),:) = [];
active_variable_matrix(active_variable_matrix(:,3) == 1 & ismember(active_variable_matrix(:,1), [3, 18, 7, 8, 11, 12]),:) = [];

[all_neurons_data, active_variable_matrix] = GetSecondStepForLightResponders(all_neurons_data, active_variable_matrix);

%%
all_active_neurons = all_neurons_data;
mean_data = nanmean(all_active_neurons, 1);
std_data = nanstd(all_active_neurons, [], 1);

% Normalize each time point in each trace
for i = 1:size(all_active_neurons,1)
    all_active_neurons(i, :) = (all_active_neurons ...
        (i, :) - mean_data) ./ std_data;
end



%%
% Remove nan rows

active_variable_matrix(sum(isnan(all_active_neurons),2) > 0,:) = [];
all_neurons_data(sum(isnan(all_neurons_data),2) > 0,:) = [];
all_active_neurons(sum(isnan(all_active_neurons),2) > 0,:) = [];

start = 20;
cols = 140;

[coef,~,~,~,explained]= pca(all_active_neurons(:,start:cols));
all_active_neurons = all_active_neurons(:,1:size(coef,1));

combinesides = 1;
[~,mean_trace_in_PC] = CalculateMeanActiveTrace(data,all_active_neurons,active_variable_matrix,coef,combinesides);



%%
% K-means and scatter plot of traces in PC space

X = mean_trace_in_PC;
PCs = [3,4];
[idx] = kmeans(X(:,PCs),3,'Replicates',500);
figure;
scatter(X(:, PCs(1)), X(:, PCs(2)), 50, idx, 'filled');
xlabel(strcat('PC',num2str(PCs(1)), {' '}, num2str(explained(PCs(1))), '%'))
ylabel(strcat('PC',num2str(PCs(2)), {' '}, num2str(explained(PCs(2))), '%'))

%%
% Plot of variance explained by PCs
figure
bar(explained)
ylabel('Variance explained')
xlabel('Principal component')

end
