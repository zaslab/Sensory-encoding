function strain_matrices = ConstructColorMatrix(unique_ids,all_acts_diffs, offstep)
strain_matrices = {};

neurons = unique(unique_ids(:,1));
for strain = 1:max(unique_ids(:,4))
    strain_matrix = [];
    strain_acts = all_acts_diffs(unique_ids(:,4) == strain);
    strain_ids = unique_ids(unique_ids(:,4) == strain,:);
    condcounter = 1;
    for cond = 1:max(unique_ids(:,2))
   
        for neuron_num = 1:length(neurons)
            neuron = neurons(neuron_num);
            strain_matrix{neuron_num, condcounter} = strain_acts{strain_ids(:,1) == neuron & strain_ids(:,2) == cond & strain_ids(:,3) == 1};
            if offstep
            strain_matrix{neuron_num, condcounter + 1} = strain_acts{strain_ids(:,1) == neuron & strain_ids(:,2) == cond & strain_ids(:,3) == offstep+1};
            end
            
        end
        condcounter = condcounter + 2;
    end
    strain_matrices{strain} = strain_matrix;
end
end