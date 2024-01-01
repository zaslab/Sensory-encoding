function confusions = GetConfusions(test_strain, condition, set_name)
letters = char(97:150);
load(strcat('train_WT_test_',test_strain, condition,...
    set_name, 'confusions.mat'))
names = who;
counter = 1;

for num = 1:length(names)
    if contains(letters, names{num})
        confusions{counter} = eval(names{num});
        counter = counter + 1;
    end
end
end