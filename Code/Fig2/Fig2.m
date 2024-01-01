function Fig2(data)
% This function plots correlation matrices for individual stimuli and the
% total mean correlation for all stimuli.

load('relevant_neurons.mat')
load('ordered_neurons')

sumsortedcorr = 0;
avrgcorrs = [];

for i = 1:7

    [avrgcorrs{i}] = ...
        CalcMeanCorr(relevant_neurons, ...
        data{i,1}(:,[40:120, 160:240],:) ,'sortAWCs','title',data{i,3});
    sumsortedcorr = sumsortedcorr + avrgcorrs{i};
    colormap(brewermap([],"Blues"))
    caxis([-0.4, 1])

end

means_sorted_corr = sumsortedcorr / 7;
means_sorted_corr = means_sorted_corr(ordered_neurons,ordered_neurons);

CGobj = clustergram(means_sorted_corr, 'RowLabels', ...
    relevant_neurons(ordered_neurons),'ColumnLabels', ...
    relevant_neurons(ordered_neurons));
close hidden

neur_corr_vals = [];
for i = 1:2:21
    neur_corr_vals(i:(i+1)) = means_sorted_corr(i,i+1);
end

[~,ind] = sort(neur_corr_vals,'descend');

mean_correlations = means_sorted_corr(ind,ind);
figure
imagesc(mean_correlations)
colormap(brewermap([],"Blues"))
f = gca;
f.YTick = 1:22;
f.XTick = 1:22;
f.YTickLabels = relevant_neurons(ordered_neurons(ind));
f.XTickLabels = relevant_neurons(ordered_neurons(ind));
title('Mean total correlation')
end

function [avrgcorrs, ind] = CalcMeanCorr(rel_neurons, current, varargin)
% Calculate correlation between neurons.
% First calculates the correlations between all neurons within
% each worm and then averages over all worms in the condition
allcors = {};
curr_title = '';
for i=1:length(varargin)
    if strcmp(varargin{i},'sortAWCs')
        sortAWC = 1;
    end

    if strcmp(varargin{i},'title')
        curr_title = varargin{i+1};
    end
end
if exist('sortAWC')
    current = sortAWCs(current);
end

for i = 1:size(current,3)

    int_neurons = rel_neurons;
    meaned = current(:,2:end,i);
    nanr = find(any(isnan(meaned),2));

    meaned(nanr,:) = [];
    int_neurons(nanr) = [];
    meaned = permute(meaned, [2,1]);
    [corr_matrix, ~] = corrcoef(meaned, 'rows', 'pairwise');

    allcors{i,1} = int_neurons;
    allcors{i,2} = corr_matrix;

end

[allcors] = FillInNans(allcors,rel_neurons);

corrformean = [];
for i = 1:size(allcors,1)
    corrformean(:,:,i) = allcors{i,2};
end

avrgcorrs = nanmean(corrformean,3);
avrgcorrs(isnan(avrgcorrs)) = 0;

figure
Z = linkage(avrgcorrs);
[~,~,ind]  = dendrogram(Z,0,'colorthreshold',0.5);
sortedcorr = avrgcorrs(ind,ind);

imagesc(sortedcorr)
title(curr_title)
a = gca;
a.YTick = 1:22;
a.XTick = 1:22;
newrons = rel_neurons;
a.YTickLabels = newrons(ind);
a.XTickLabels = newrons(ind);
c = colorbar('location', 'eastoutside');
set(gca,'XTickLabelRotation',45)
colormap jet;

end

%%
function [allcors] = FillInNans(allcors,int_neurons)

% This function adds nan rows and columns to replace missing neurons in the
% correlation matrices
for i = 1:size(allcors,1)

    missN =find(ismember(int_neurons,allcors{i,1})==0);

    if isempty(missN) == 0

        NmissN = length(missN);

        for j = 1:NmissN

            allcors{i,2} = [allcors{i,2}(1:missN(j)-1,:); ...
                NaN(1,size(allcors{i,2},1));...
                allcors{i,2}(missN(j):end,:)];

            allcors{i,2} = cat(2,allcors{i,2}(:,1:missN(j)-1),...
                NaN(size(allcors{i,2},1),1),...
                allcors{i,2}(:,missN(j):end));

        end
    end
end
end
