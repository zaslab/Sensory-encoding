function [bothsteps_matrix, bothstepids] = MakeWormFeatureVectors(strain_matrices,offstep)
asf
strain_matrix = [];
id_matrix = [];
for strain = 1:length(strain_matrices)
%     curr_strain = RemoveCrossReads(strain_matrices{strain}, strain);

    curr_strain = strain_matrices{strain};
    neuron_means = cellfun(@nanmean, curr_strain);

    for cond = 1:size(neuron_means,2)
        curr_cond = curr_strain(:,cond);
        cond_means = neuron_means(:,cond);
        for worm = 1:size(curr_cond{1},1)
            curr_worm = cellfun(@(v)v(worm),curr_cond)';
            curr_worm(isnan(curr_worm)) = cond_means(isnan(curr_worm));
            strain_matrix = cat(1,strain_matrix, curr_worm);
            id_matrix = cat(1,id_matrix, [worm, cond,strain]);
        end
    end
end
bothstepids = [];
bothsteps_matrix = [];
for i = 1:max(id_matrix(:,2))
    
    for strain = 1:3
        if strain ~= 1 & i > 7
            continue
        end
        onsteps = strain_matrix(id_matrix(:,2) == i & id_matrix(:,3) == strain,:);
        if offstep
        offsteps = strain_matrix(id_matrix(:,2) == i + 1 & id_matrix(:,3) == strain , :);
        
        bothsteps = cat(2,onsteps,offsteps);
        else
            bothsteps = onsteps;
        end
        ids = [i * ones(size(bothsteps,1),1), strain * ones(size(bothsteps,1),1)];
        bothstepids = cat(1,bothstepids,ids);
        bothsteps_matrix = cat(1,bothsteps_matrix,bothsteps);
%         method = 'znorm';
        unnormed_bothsteps_matrix = bothsteps_matrix;
%         bothsteps_matrix = NormalizeData(unnormed_bothsteps_matrix, method);
    end
end
bothstepids(:,3) = bothstepids(:,1) * 10 + bothstepids(:,2);

end