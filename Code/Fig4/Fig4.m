function Fig4(WTdata, unc31data, unc13data)
load('relevant_neurons.mat')
load('ordered_neurons.mat')

data = CleanSensoryData(WTdata, unc31data, unc13data);

all_neurons_data = [];
variable_matrix = [];
for strain = 1:length(data)
    [curr_neurons_data, curr_variable_matrix] = prepare_PCA_data(data{strain}, 'step_benchmark',strain);
    all_neurons_data = cat(1, all_neurons_data, curr_neurons_data);
    variable_matrix= cat(1, variable_matrix, curr_variable_matrix(:,1:4));
end

for u = 1:size(all_neurons_data,1)
    all_neurons_data(u,:) = smooth(all_neurons_data(u,:));
end

all_neurons_data = Reallign(all_neurons_data, variable_matrix);

[all_neurons_data, active_variable_matrix] = GetSecondStepForLightResponders(all_neurons_data, variable_matrix);
[~, unique_ids] = calculate_step_means(active_variable_matrix, all_neurons_data, 0);

% Currently: WT are compared to rhodamine in a two-group ttest, mutants are
% compared to WTs. Correction performed on each strain separately
NeuronsUsed = {'ASI', 'AWCoff', 'AWCon', 'ASJ', 'ADL', 'ASK','ASH','ADF','AWA','AWB','ASG','ASER','ASEL'};

for u = 1:size(all_neurons_data,1)
    all_neurons_data(u,:) = smooth(all_neurons_data(u,:));
end
%%
all_act_diffs = calculate_activations(active_variable_matrix, all_neurons_data, 0,...
    'all_window', 20:50, 'ASH_window', 26:40);

offstep=1;
strain_matrices = ConstructColorMatrix(unique_ids,all_act_diffs,offstep);


for strain = 1:3
    strain_matrices{strain} = CompleteStrainMatrix(strain_matrices{strain});
end

colorcelldata = ConstructColorAndPvalMatrices(strain_matrices);
num_of_neurons = size(colorcelldata(1).activations,1);
% Fill out the missing field for convenience
colorcelldata(1).MutToWTpvalues = ones(num_of_neurons,14);
corrected = 'fdr';
for i = 1:3

    ctrlpvals = colorcelldata(i).ctrlpvalues;
    MutToWTpvalues = colorcelldata(i).MutToWTpvalues;
    corrected_ps_data(i).activations= colorcelldata(i).activations;
    corrected_ps_data(i).ctrlpvalues = ctrlpvals;
    corrected_ps_data(i).MutToWTpvalues = MutToWTpvalues;

end
%%
WTstd = colorcelldata(1).stds;
WTacts = colorcelldata(1).activations;

% Remove significance if WT and mut are less than WTstd apart 

for i = 2:3
    std_diff = ones(num_of_neurons,14);
    mutacts = colorcelldata(i).activations;
    std_diff(abs(WTacts - mutacts) < WTstd) = 0;
    corrected_ps_data(i).bigger_than_std = std_diff;
end
corrected_ps_data(1).bigger_than_std = zeros(num_of_neurons,14);
%%
color_matrix = {};
for type = 1:4
    color_matrix{type} = CombineStrainMatrices(corrected_ps_data,type);
end

act_matrix = color_matrix{1};
ctrlpvals_matrix = color_matrix{2};

MutToWTpvalues_matrix = color_matrix{3};
bigger_than_std_matrix = color_matrix{4};

ctrlpvals_matrix  = MultCompCorrection(ctrlpvals_matrix ,corrected);
MutToWTpvalues_matrix = MultCompCorrection(MutToWTpvalues_matrix,corrected);

pthresh = 0.05;
amp_thresh = 0.1;
MutToWTpvalues_matrix(bigger_than_std_matrix == 0) = 1;
orig_mut_pvals = MutToWTpvalues_matrix;
MutToWTpvalues_matrix(:,1:3:end) = 0;
color_matrix_p_accounter = act_matrix;
ctrlpvals_matrix(color_matrix_p_accounter == 0) = 1;


DrawColorCells(color_matrix_p_accounter, 'labels', NeuronsUsed);
AddSignificance(ctrlpvals_matrix,MutToWTpvalues_matrix,color_matrix_p_accounter,pthresh,amp_thresh);

end

function AddSignificance(ctrlpvals_matrix,MutToWTpvalues_matrix,color_matrix_p_accounter,pthresh,amp_thresh)
for row = 1:size(ctrlpvals_matrix,1)
    for col = 1:size(ctrlpvals_matrix,2)
        strain = mod(col,3);
        if strain ~= 1
            if strain == 2
                WT_modifier = 1;
            else
                WT_modifier = 2;
            end
            if MutToWTpvalues_matrix(row, col) < pthresh &...
                    ~(ctrlpvals_matrix(row, col) > pthresh &...
                    ctrlpvals_matrix(row, col - WT_modifier) > pthresh)
                text(col-0.2,row+0.1,'*', 'FontSize', 30, 'Color', 'k')
            
            end
        else
            if ctrlpvals_matrix(row, col) < pthresh & abs(color_matrix_p_accounter(row, col)) > amp_thresh
                text(col-0.2,row+0.1,'*', 'FontSize', 30, 'Color', 'm')
            end
        end
    end
end
end
%%
function all_data = ConstructColorAndPvalMatrices(strain_matrices)

load("neuron_pairs.mat")
num_of_neurons = size(neuron_pairs,1);
swpvals_matrix = [];
for strain = 1:3

    [WTref,strain_ref] = GetStrainControl(strain_matrices,strain);
    current_strain = strain_matrices{strain};

    current_color_matrix = nan(num_of_neurons,12);
    current_std_matrix = nan(num_of_neurons,12);

    for cond = 1:2:12
        ref_acts = [];
        ref_pval = [];

        for neuron = 1:size(neuron_pairs,1)
            current_neuron_act_steps = [];

            current_neuron_steps = GetCurrentSteps(current_strain,cond,neuron_pairs,neuron);
            ref_steps = GetRefSteps(strain_ref,neuron_pairs,neuron);


            current_neuron_act_steps(1,:) = current_neuron_steps(1,:);% - nanmean(ref_steps(1,:));
            current_neuron_act_steps(2,:) = current_neuron_steps(2,:);% - nanmean(ref_steps(2,:));


            current_color_matrix(neuron,cond) = nanmean(current_neuron_act_steps(1,:));
            current_color_matrix(neuron,cond + 1) = nanmean(current_neuron_act_steps(2,:));

            current_std_matrix(neuron,cond) = nanstd(current_neuron_act_steps(1,:))/sqrt(sum(~isnan(current_neuron_act_steps(1,:))));
            current_std_matrix(neuron,cond + 1) = nanstd(current_neuron_act_steps(2,:))/sqrt(sum(~isnan(current_neuron_act_steps(1,:))));

            ref_pval(neuron,1:2) =  OneSampleTtest(ref_steps);

            [~,WT_control] = GetStrainControl(strain_matrices,1);
            WT_control_steps = GetRefSteps(WT_control,neuron_pairs,neuron);
            current_WT_steps = GetCurrentSteps(WTref,cond,neuron_pairs,neuron);

            if strain == 1
                %
                control_pvals = TwoSampleTtest(ref_steps, current_neuron_act_steps);
                [~, swpvals(1)] = swtest(current_neuron_act_steps(1,:));
                [~, swpvals(2)] = swtest(current_neuron_act_steps(2,:));
            else
                % Compare mutants to both rhodamine control and WT
                [control_pvals,WT_pvals] = GetMutantPvals(current_neuron_act_steps,current_WT_steps);
                mutant_to_WT_cond_pvals(neuron,cond:(cond + 1)) = WT_pvals;

                mutant_ctrl_pvals(neuron,:)  = TwoSampleTtest(ref_steps, WT_control_steps);


            end

            control_pval_matrix(neuron,cond:(cond + 1)) = control_pvals;
            swpvals_matrix(neuron,cond:(cond + 1)) = swpvals;
            ref_acts(neuron,[1 2]) = nanmean(ref_steps,2);
            ref_std(neuron,[1 2]) = nanstd(ref_steps,0,2)./sqrt(sum(~isnan(ref_steps),2));

        end
    end

    strain_color_matrix = cat(2,ref_acts,current_color_matrix);
    strain_pval_matrix =  cat(2,ref_pval,control_pval_matrix);
    strain_std_matrix =   cat(2,ref_std,current_std_matrix);

    all_data(strain).activations = strain_color_matrix;
    all_data(strain).ctrlpvalues = strain_pval_matrix;
    all_data(strain).stds = strain_std_matrix;

    if strain > 1
        mutant_to_WT_pvals = cat(2,mutant_ctrl_pvals,mutant_to_WT_cond_pvals);
        all_data(strain).MutToWTpvalues = mutant_to_WT_pvals;
    end


end
end

%%

function color_matrix = CombineStrainMatrices(data,type)

types = {'activations','ctrlpvalues','MutToWTpvalues','bigger_than_std'};
current_type = types{type};
color_matrix = nan(size(data(1).activations,1),42);

for strain = 1:3

    current_matrix = data(strain).(current_type);

    for i = 1:14
        color_matrix(:,(i - 1) * 3 +  strain) = current_matrix(:,i);
    end


end
end
%%
function corrected_ps = MultCompCorrection(pvals, corrected)

FDR_pvals = mafdr(pvals(:),'BHFDR', true);
FDR_pvals = reshape(FDR_pvals',size(pvals,1),size(pvals,2));

[p,~,~,~] = fwer_bonf(pvals, 0.05, false);
bonf_pvals = reshape(p,size(pvals,1),size(pvals,2));

if strcmp(corrected, 'fdr')
    corrected_ps = FDR_pvals;
elseif strcmp(corrected,'bonf')
    corrected_ps = bonf_pvals;
end

end
%%
function [WT, strain_ref] = GetStrainControl(strain_matrices,strain)
% Returns the rhodamine control responses for the requested strain, as well
% as the entire WT response matrix
strain_ref = strain_matrices{strain}(:,[13,14]);
WT = strain_matrices{1};
end


%%
function ref_steps = GetRefSteps(strain_ref,neuron_pairs,neuron)
% Extracts activation values of the rhodamine control of the requested
% neuron and strain. Combine left and right sides for all except for ASE
% and AWC

% OUTPUT -2 row matrix. first row - onstep activations, scnd row - off step

neuron1 = neuron_pairs(neuron,1);
neuron2 = neuron_pairs(neuron,2);
%Check if AWC or ASE.

if isnan(neuron_pairs(neuron,2))
    ref_on_step = strain_ref{neuron1,1};
    ref_off_step = strain_ref{neuron1,2};
else
    ref_on_step = cat(1,strain_ref{neuron1,1}, strain_ref{neuron2,1});
    ref_off_step = cat(1,strain_ref{neuron1,2}, strain_ref{neuron2,2});
end

ref_steps = [ref_on_step, ref_off_step]';
end

%%
function step_pvals = TwoSampleTtest(current_neuron_act_steps, ref_adjusted_steps)

[~, pval_onstep] = ttest2(current_neuron_act_steps(1,:),ref_adjusted_steps(1,:));
[~, pval_offstep] = ttest2(current_neuron_act_steps(2,:),ref_adjusted_steps(2,:));

step_pvals = [pval_onstep,pval_offstep];

end

%%
function [control_pvals,WT_pvals] = GetMutantPvals(neuron_act_steps,WT_steps)

% For mutant test if they are different from WT and also if they are
% different from 0

control_pvals =  OneSampleTtest(neuron_act_steps);

WTref_adjusted_steps(1,:) = WT_steps(1,:);
WTref_adjusted_steps(2,:) = WT_steps(2,:);

WT_pvals  = TwoSampleTtest(neuron_act_steps, WTref_adjusted_steps);

end

%%
function current_neuron_steps = GetCurrentSteps(current_strain,cond,neuron_pairs,neuron)
% Extracts the activation values of the requested neuron. Combine left and
% right sides for all except for ASE and AWC

% OUTPUT -2 row matrix. first row - onstep activations, scnd row - off step

neuron1 = neuron_pairs(neuron,1);
neuron2 = neuron_pairs(neuron,2);

%Check if AWC or ASE.
if isnan(neuron_pairs(neuron,2))
    current_neuron_on_step = current_strain{neuron1, cond};
    current_neuron_off_step = current_strain{neuron1, cond + 1};
else
    current_neuron_on_step = cat(1,current_strain{neuron1, cond},current_strain{neuron2, cond});
    current_neuron_off_step = cat(1,current_strain{neuron1, cond + 1},current_strain{neuron2, cond + 1});
end

current_neuron_steps = [current_neuron_on_step,current_neuron_off_step]';
end
%%
function step_pvals =  OneSampleTtest(steps)

[~,on_pval] = ttest(steps(1,:));
[~,off_pval] = ttest(steps(2,:));

step_pvals = [on_pval,off_pval];
end

