function Fig5()
conditions = {'', '_volatility_', '_valence_',''};
cond_names = {{'IAA' 'DA' 'Gly' 'NaCl' 'Quin' 'SDS' 'Ctrl'}, {'Soluble', 'Volatile'}, {'Attractive','Aversive'}};
test_strain = {'WT','31','13'};
set_names = {'both_steps_all_'};

for strain = 1:3
    for cond = 1:3

        set_name = [conditions{cond}, set_names{1}];
        confusions = GetConfusions(test_strain{strain}, '', set_name);
        confusions_cat = cat(3,confusions{:})./31;
        PlotConfusions(confusions_cat , cond_names{cond}, test_strain{strain})

    end
end
end

function [weighted_F,  all_metrics] = PlotConfusions(confusions, cond_names,strain)
load('rocket.mat')
neurons  = {'AWCR','ASJ','ASH','ASG','ADF','ASI','ADL','ASK','ASER','AWA','AWB','AWCL','ASEL'};
all_metrics = [];

for j = 0:4:12
    figure
    for counter = 1:4

        i= counter + j;

        if i > 13
            continue
        end

        neuron_ids = i;

        subplot(2,2,counter)
        imagesc(confusions(:,:,i))
        title(neurons{neuron_ids(1)})
        colormap(rocket)
        caxis([0,1])
        set(gca, 'XTick', 1:length(cond_names))
        set(gca, 'YTick', 1:length(cond_names))
        set(gca, 'XTickLabel', cond_names)
        set(gca, 'YTickLabel', cond_names)
        xtickangle(45)

    end
    sgtitle(strain)
end

end

