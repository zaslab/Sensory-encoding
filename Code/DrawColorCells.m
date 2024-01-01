function DrawColorCells(Data, varargin)
for i = 1:2:length(varargin)
    if strcmp('labels', varargin{i})
        labels = varargin{i+1};
    else
        labels = {'ASI', 'AWCL','AWCR','ASJ','ASH','AWA','AWB','ASER','ASEL'};
    end
end
figure
imagesc(Data)
set(gca,'XTick', 2:3:42, 'XTickLabel',{'rho On', 'rho Off', 'IAA On', 'IAA Off', 'DA On', 'DA Off', 'Gly On', 'Gly Off',...
    'NaCl On', 'NaCl Off','Qui On', 'Qui Off','SDS On', 'SDS Off'})
set(gca,'YTick', 1:length(labels), 'YTickLabel',labels)


[rows, columns, ~] = size(Data);
hold on;
for row = 0.5 : 1 : rows + 5
    line([0, columns + 5], [row, row], 'Color', 'k', 'LineWidth', 1);
end
for col = 0.5 : 3 : (columns +3)
    line([col, col], [0, rows + 5], 'Color', 'k', 'LineWidth', 1);
end
MakeColorCellsGreenBlue
end


function MakeColorCellsGreenBlue

blue = [0.88, 0.67, 0.004];
yellowstart = [1 1 1];
yellow =  [0 0 1];
blueend =[1, 1, 1];

Length = 175;
colors_p = [linspace(yellow(1),blue(1),Length)', ...
    linspace(yellow(2),blue(2),Length)', linspace(yellow(3),blue(3),Length)'];
greenend = colors_p(end,:);
Length = 123;
Start = greenend;
End = yellowstart;
yellowmap= [linspace(End(1),Start(1),Length)',...
    linspace(End(2),Start(2),Length)', linspace(End(3),Start(3),Length)'];


bluestart = colors_p(1,:);
Length = 52;
Start = blueend;
End = bluestart;
bluemap = [linspace(End(1),Start(1),Length)',...
    linspace(End(2),Start(2),Length)', linspace(End(3),Start(3),Length)'];
allcolormap = cat(1, bluemap,yellowmap);
colormap(allcolormap)
caxis([-0.5, 1.2])
end
