function CompletedMatrix = CompleteStrainMatrix(strain_matrix)
CompletedMatrix = {};
for cond = 1:size(strain_matrix,2)
    cond_matrix = [];
    for neuron = 1:size(strain_matrix,1)
        cond_matrix(neuron,:) =strain_matrix{neuron, cond};
    end
        completed_cond_matrix = CompleteMatrix(cond_matrix);
        for completed_neuron = 1:size(strain_matrix,1)
            CompletedMatrix{completed_neuron,cond} = completed_cond_matrix(completed_neuron,:)';
        end
end
end