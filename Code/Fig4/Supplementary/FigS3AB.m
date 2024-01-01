%Plot ASJs

function FigS3AB(WTdata, unc31data)
onstep = 60:180;
offstep = 180:300;

cond1 = WTdata{1,1};
cond2 = unc31data{1,1};

[WTASJON,WTASJOFF,WTASJactsON,WTASJactsOFF, WTASJacts, WTASJsIAA]...
                                                 = GetASJActs(WTdata{1,1});

[unc31ASJON,unc31ASJOFF,unc31ASJactsON,unc31ASJactsOFF,unc31ASJacts, unc31ASJsIAA]...
                                                = GetASJActs(unc31data{1,1});
plotASJheatmap(WTASJsIAA)
title('WT')

plotASJheatmap(unc31ASJsIAA)
title('unc-31')

meanacts = [nanmean(WTASJacts),nanmean(unc31ASJacts)];
actste = [std(WTASJacts)/sqrt(length(WTASJacts)),...
    std(unc31ASJacts)/sqrt(length(unc31ASJacts))];

meanacts = [WTASJON, WTASJOFF, unc31ASJON, unc31ASJOFF];

actste = [std(WTASJactsOFF)/sqrt(length(WTASJactsON)),...
    std(WTASJactsOFF)/sqrt(length(WTASJactsOFF)),...
    std(unc31ASJactsON)/sqrt(length(unc31ASJactsON)),...
    std(unc31ASJactsOFF)/sqrt(length(unc31ASJactsOFF))];

[~, PvalOn] = ttest2(WTASJactsON, WTASJactsOFF);
[~, PvalOff] = ttest2(unc31ASJactsON, unc31ASJactsOFF);

plotASJbar(meanacts,actste)

hold on;
UnivarScatter(padcat(WTASJactsON,WTASJactsOFF,unc31ASJactsON,unc31ASJactsOFF), ...
    'Whiskers', 'lines');

end
 
function [ASJON,ASJOFF,actsON, actsOFF, acts,ASJs] = GetASJActs(data)

onstep = 60:180;
offstep = 180:300;

ASJR = permute(data(2,2:301,:),[3,2,1]);
ASJL = permute(data(20,2:301,:),[3,2,1]);
ASJs = cat(1,ASJL,ASJR);

mean_data = nanmean(ASJs,2);
ASJs(isnan(mean_data),:) = [];
mean_data(isnan(mean_data)) = [];

[~, idx] = sort(mean_data,'descend');
ASJs = ASJs(idx,:);

acts1 =  mean(ASJs(:,onstep),2) - mean(ASJs(:,1:60),2);
acts2 =  mean(ASJs(:,offstep),2) - mean(ASJs(:,onstep),2);

actsON = nanmean([ASJs(:,onstep)],2);
actsOFF = nanmean([ASJs(:,offstep)],2);

ASJON = mean(actsON);
ASJOFF = mean(actsOFF);
acts = [acts1',acts2'];

end

function plotASJheatmap(data)
figure

imagesc(data)
colormap jet
caxis([1, 2.2])
ax = gca;
ax.YTick = 0:5:50;
ax.XTick = 60:60:299;
ax.XTickLabel = 30:30:240;
ax.FontSize = 12;
ylabel('Worms')
end

function plotASJbar(meanacts, actste)

figure
hold on

h = bar(meanacts,'BarWidth',0.4, 'FaceColor', 'flat');
nbars = size(meanacts, 2);

% Calculating the width for each bar group
for i = 1:nbars
    x = i;
    errorbar(x, meanacts(:,i), actste(:,i), 'k.');
end

ylim([1 max(meanacts(:)) + max(actste(:)) ])
set(gca,'xtick',[])
legend({'ON step', 'OFF step'}, 'FontSize', 15)
ylim([0.95, 1.53])
legend off

end