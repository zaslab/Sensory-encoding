function cleaned_data = CleanSensoryData(WT, unc31, unc13)

% Sort AWCs by strength
for i = 1:7
    if i == 4
        step = 1;
    else
        step = 2;
    end
    WT(i,:) = sortByStrength(WT(i,:), 1,{'AWC'},step);
     unc13(i,:) = sortByStrength(unc13(i,:), 1,{'AWC'},step);
    unc31(i,:) = sortByStrength(unc31(i,:), 1,{'AWC'},step);
end


   cleaned_data = {WT, unc31, unc13};
end


