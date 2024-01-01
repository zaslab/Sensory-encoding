function Fig1(data)
load('relevant_neurons.mat')
order = [13,5,12,7,21,9,19,4,18,3,15,6,20,2,11,8,17,10,14,22,16,1];

for cond = 1:7
matrix = sortAWCs(data{cond,1});
step_locations = data{cond,4};
for step = 1:2
    figure
    current_data = zeros(22,120,size(matrix,3));
    for neuron = 1:22
        for worm = 1:size(matrix,3)
            step_loc = step_locations(neuron, step, worm);
            if step_loc > 300
                step_locations(neuron, step, worm) = nan;
                step_loc = nanmean(step_locations(:,step,worm));
            end
            current_data(neuron, :, worm) = matrix(neuron, (step_loc-20):(step_loc+99), worm);
        end
    end
    imagesc(nanmean(current_data(order,:,:),3))
    hold on
    plot([20 20], [0 25], 'w', 'LineWidth',3)
    yticks(1:22)
    yticklabels(relevant_neurons(order))
    colormap jet
   
    caxis([1, 2.25])
    title(data{cond,3})
end
end

           
end