%%
function FigS3EF(WTdata, unc31data, unc13data)
load('relevant_neurons.mat')
load('ordered_neurons.mat')

cond = 3;
neurons = [3,18];


for step = 1:2
    if step == 1
        window = 2:30;
    else
        window = 2:100;
    end
    
    alignedWT = AlignToPeak(WTdata, step, neurons, cond);
    alignedunc13= AlignToPeak(unc13data, step, neurons, cond);
    alignedunc31= AlignToPeak(unc31data, step, neurons, cond);

    meanWTtrace = smooth(nanmean(alignedWT,1))';
    meanWTtrace = meanWTtrace - min(meanWTtrace(window));
    errorsWT = nanstd(alignedWT,1)./sqrt(size(alignedWT,1));
    meanunc31trace = nanmean(alignedunc31,1);
    meanunc31trace  = meanunc31trace - min(meanunc31trace(window));
    errorsunc31 = nanstd(alignedunc31,1)./sqrt(size(alignedunc31,1));
    meanunc13trace = nanmean(alignedunc13,1);
    meanunc13trace  = meanunc13trace - min(meanunc13trace(window));
    errorsunc13 = nanstd(alignedunc13,1)./sqrt(size(alignedunc13,1));

    figure
    hold on
    x = 1:length(meanWTtrace);
    h1 = fill([x fliplr(x)],...
        [meanWTtrace-errorsWT fliplr(meanWTtrace+errorsWT)]...
                                        ,'b', 'linestyle', 'none');
    alpha(.2)
    h2 =  plot(x,meanWTtrace,'b', 'LineWidth',3);

    h3 = fill([x fliplr(x)],...
        [meanunc31trace-errorsunc31 fliplr(meanunc31trace+errorsunc31)],...
                                             'r', 'linestyle', 'none');
    alpha(.2)
    h4 = plot(x,meanunc31trace,'r', 'LineWidth',3);

    h5 = fill([x fliplr(x)],...
        [meanunc13trace-errorsunc13 fliplr(meanunc13trace+errorsunc13)],...
                                             'k', 'linestyle', 'none');
    alpha(.2)
    h6 = plot(x,meanunc13trace,'k', 'LineWidth',3);

    ax = gca;
    ax.XTick = 25:20:480;
    ax.XTickLabels = {0 10 20 30 40};
    xlim([15 85])
    legend([h2 h4 h6],'WT','unc-31','unc-13')

end

end

%%
function aligned= AlignToPeak(datastruct, step, neurons, cond)

data = permute(cat(3,datastruct{cond ,1}(neurons(1),:,:),...
    datastruct{cond ,1}(neurons(2),:,:)),[3,2,1]);

steps = permute(cat(3,datastruct{cond ,4}(neurons(1),step,:),...
    datastruct{cond ,4}(neurons(2),step,:)),[3,2,1]);

aligned = [];
ind = [];
for i = 1:size(data,1)
    window = (steps(i) - 18):(steps(i) + 15);
    [~, ind(i)] = max(data(i,window));
    ind(i) = ind(i) + window(1) - 1;
    aligned(i,:) = data(i, (ind(i) - 30):(ind(i) + 80));
    aligned(i,:) = aligned(i,:) - mean(aligned(i,1:10));

end
end



