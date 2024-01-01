function Fig6EF()
accuracy = CalculateAccuracy();

flag = 'Lin';
names = GetNames(flag);
for i = 1:size(names,1); name_cell{i} = strrep(names(i,:),'synth_Lin_','');end
for i = 1:size(names,1); name_cell{i} = strrep(name_cell{i},'_','');end
for i = 1:size(names,1); name_cell{i} = strrep(name_cell{i},'_','');end
for i = 1:size(names,1); name_cell{i} = strrep(name_cell{i},'.mat','');end;
for i = 1:size(names,1); name_cell{i} = strrep(name_cell{i},' ','');end
name_cell = name_cell';
%%
use_names = [1:9];

colors = colormap(jet(66));
close
use_colors = [colors(15,:);
    [0.2 0.9 0.2];
    colors(62,:)];
use_colors = repmat(use_colors, 3,1);

figure; UnivarScatter(accuracy(use_names,:)',  'MarkerFaceColor', use_colors, 'Whiskers','lines')
set(gca, 'XTickLabel', name_cell(use_names))
set(gca, 'Ytick', 0:0.2:1)
ylim([0 1])
[pvals, bonf] = GetLinPvals(accuracy(use_names,:)');
%%

use_names = [10:12];

figure; UnivarScatter(accuracy(use_names,:)', 'PointSize', 60, 'MarkerFaceColor', use_colors, 'Whiskers','lines')
set(gca, 'XTickLabel', name_cell(use_names))
set(gca, 'Ytick', 0:0.2:1)
ylim([0.1 0.62])
xlim([0 4])
end
%%
function accuracy = CalculateAccuracy()

flag = 'synth_Lin';
names = GetNames(flag);
accuracy = [];

for run = 1:size(names,1)

    clearvars -except run accuracy
    flag = 'synth_Lin';
    names = GetNames(flag);
    load(strrep(names(run,:),' ',''))
    clear names
    clear flag
    names = who;
    names(strcmp(names, 'run')) = [];
    names(strcmp(names, 'accuracy')) = [];
    counter = 1;

    for num = 1:length(names)
        confusions{counter} = eval(names{counter});
        counter = counter + 1;

    end

    confusionMatrix = confusions{1};
    conf_repeats  = [];

    for rep = 1:length(confusionMatrix)
        current_matrix = confusionMatrix{rep};
        conf_repeats = cat(3, conf_repeats, current_matrix );
        accuracy(run,rep) = sum(diag(current_matrix ))/sum(current_matrix(:));
    end
end
end
%%

function names = GetNames(flag)
all_names = ls;
counter  = 1;
for i = 1:size(all_names,1)

    if contains(all_names(i,:), flag)

        names(counter,:) = all_names(i,:);
        counter = counter + 1;

    end
end
end
%%

function [pvals, FDR] = GetLinPvals(data)
counter = 1;
pvals = [];

for i = 1:8
    if ismember(i,[3,6])
        continue
    end

    for j = i:(i + 2)
        if j == i
            continue
        end
        if ismember(i, [2,5]) & ismember(j,[4,6])
            continue
        end
        if j == 10
            continue
        end

        pvals(counter,1:2) = [i,j];
        [~, pvals(counter,3)] = ttest2(data(:,i), data(:,j));
        counter = counter + 1;
    end
end

pvals(end+1,1:2) = [1,4];
[~, pvals(end,3)] = ttest2(data(:,1), data(:,4));
pvals(end+1,1:2) = [1,7];
[~, pvals(end,3)] = ttest2(data(:,1), data(:,7));
pvals(end+1,1:2) = [4,7];
[~, pvals(end,3)] = ttest2(data(:,4), data(:,7));
FDR = mafdr(pvals(:,3)','BHFDR',1);

end