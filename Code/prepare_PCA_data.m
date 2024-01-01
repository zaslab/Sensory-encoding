function [data_for_pca, variable_matrix] = prepare_PCA_data(data, activation_start,strain)
% Allign activations to either the step (earliest response), activation, or
% step benchmark determined by fastest responding neurons
% of each neuron - 'step', 'activation', 'step_benchmark'

sides = [1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 1];
if strcmp(activation_start, 'activation')
    end_idx = 121;
else
    end_idx = 180;
end
variable_matrix = [];
data_for_pca = [];

counter = 1;

current_data = data;
mat_sizes = cellfun(@(x) size(x,3) , current_data(:,1), 'UniformOutput', false);
step_window = [40, 140];
wormnum = sum(cell2mat(mat_sizes));
current_variable_matrix = nan(wormnum*22*2,6);
variable_matrix = cat(1, variable_matrix, current_variable_matrix);
current_all_neurons_data = nan(wormnum*22*2, end_idx);
data_for_pca = cat(1, data_for_pca, current_all_neurons_data);

for cond = 1:7
    current_cond = current_data{cond,1};
    current_steps = current_data{cond, 4};
    for worm = 1:size(current_cond,3)
        current_worm = current_cond(:,:,worm);
        for step = 1:3
            if strcmp(activation_start, 'step')
                earliest_response = repmat(min(current_steps(:,step,worm)),[22,1]);


            elseif strcmp(activation_start, 'activation')
                earliest_response = current_steps(:,step,worm);

            elseif strcmp(activation_start, 'step_benchmark')
                earliest_response = BenchMarkStepStarts(current_steps(:,step, worm),cond, step);

            end
            for neuron = 1:22
                % clean miss-marked step starts to avrg of all other
                % neurons

                [earliest_response, current_worm(neuron,:)] = CheckStepValidity(earliest_response,neuron, current_worm(neuron,:),step_window(2));
                current_step(neuron,:) = current_worm(neuron, ...
                    (earliest_response(neuron) - step_window(1)):(earliest_response(neuron) + step_window(2)));
       
                if strcmp(activation_start, 'activation')
                    [~, idx] = max(current_step(neuron, 20:80));
                else
                    idx = 1;
                end


                current_step_sub(neuron,:) = current_step(neuron,idx:(idx + end_idx - 1));
            end
            data_for_pca(counter:(counter + 21),:) = current_step_sub;

            variable_matrix(counter:(counter + 21),:) = cat(2,(1:22)', repmat([cond, step],22,1), repmat(strain,22,1),repmat(worm,22,1), sides');

            counter = counter + 22;
        end
    end
end


end

function fixed_step_starts = clean_step_starts(all_step_starts, neuron_num)
fixed_step_starts = all_step_starts;
all_step_starts(neuron_num) = [];
fixed_step_starts(neuron_num) = round(nanmean(all_step_starts));
end

function step_starts = BenchMarkStepStarts(step_starts,cond, step)
if step == 3
    step = 1;
end
% IAA
if cond == 1
    if step == 1
        on_step_benchmark = [10,17,14,22];
    elseif step == 2
        off_step_benchmark = [1];
    end
end

% DA
if cond == 2
    if step == 1
        on_step_benchmark = [10,17,16,1];
    elseif step == 2
        off_step_benchmark = [16,1,14,22];
    end
end

% Gly
if cond == 3
    if step == 1
        on_step_benchmark = [3,18,16,1];
    elseif step == 2
        off_step_benchmark = [16,1,6,15];
    end
end

% NaCl
if cond == 4
    if step == 1
        on_step_benchmark = [21,22,14];
    elseif step == 2
        off_step_benchmark = [9,3,18];
    end
end

% Quin
if cond == 5
    if step == 1
        on_step_benchmark = [21,2,20];
    elseif step == 2
        off_step_benchmark = [14,22,9];
    end
end
% SDS
if cond == 6
    if step == 1
        on_step_benchmark = [21];
    elseif step == 2
        off_step_benchmark = [9,14,22];
    end
end

% DA + NaCl
if cond == 8
    if step == 1
        on_step_benchmark = [21,10,17,14,22];
    elseif step == 2
        off_step_benchmark = [9,16,1];
    end
end

% rho
if cond == 7
    if step == 1
        on_step_benchmark = [9];
    elseif step == 2
        off_step_benchmark = [14,22,1];
    end
end

% IAA + DA
if cond == 9
    if step == 1
        on_step_benchmark = [10,17,14,22];
    elseif step == 2
        off_step_benchmark = [16,1];
    end
end

% Quin + NaCl
if cond == 10
    if step == 1
        on_step_benchmark = [21,14,22];
    elseif step == 2
        off_step_benchmark = [9,6,15,1,16];
    end
end

% DA + SDS
if cond == 11
    if step == 1
        on_step_benchmark = [10,17,21];
    elseif step == 2
        off_step_benchmark = [14,22,9];
    end
end

if step == 1
    step_marks = on_step_benchmark;
elseif step == 2
    step_marks = off_step_benchmark;
end

actual_step = floor(nanmean(step_starts(step_marks)));

step_starts(:) = actual_step;
end

function [earliest_response, current_worm] = CheckStepValidity(earliest_response,neuron, current_worm,step_window)
if earliest_response(neuron) > (479 - step_window)
    earliest_response = clean_step_starts(earliest_response, neuron);

end

if earliest_response(neuron) > (479 - step_window) | isnan(earliest_response(1))
    earliest_response(:) = 100;
    current_worm = nan;
end

end

