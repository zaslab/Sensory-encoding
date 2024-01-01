function [sortedData] = sortByStrength(alldata,cond,neurons,step)
% This function sorts a neuron pair by strength. By default the R neuron 
% becomes the stronger of the pair, and the L becomes the weaker.
% If one of the neurons doesn't exist replaces both with nan
allneurons = {'AWCR','ASJR','ASHR','ASGR','ADFR','ASIR','ADLR','ASKR','ASER','AWAR',...
    'ASKL','ADLL','ADFL','AWBL','ASIL','AWCL','AWAL','ASHL','ASGL','ASJL','ASEL','AWBR'};

if exist('step') == 1
    if length(step) == 1
        window = [(60+120*(step-1)-20):(60+120*(step-1)+40)];
    else
        window = step;
    end
else
    window = 1:size(alldata,2);
end
if iscell(alldata) && size(alldata,1) > 1
    steplocs = alldata{cond,4};
    data = alldata{cond,1};
elseif iscell(alldata)
    steplocs = alldata{1,4};
    data = alldata{1,1};
else
    data = alldata;
end


for j = 1:length(neurons)
    currneur = neurons{j};
    Rind = strfind(allneurons, strcat(currneur,'R'));
    Rside = find(not(cellfun('isempty',Rind)));
    
    Lind = strfind(allneurons, strcat(currneur,'L'));
    Lside = find(not(cellfun('isempty',Lind)));
    
    
for i = 1:size(data,3)
    if isempty(Lside) || isempty(Rside)
        data(Rside,:,i) = NaN;
        data(Lside,:,i) = NaN;
        
    elseif sum(isnan(data(Rside,:,i))) > 10 || sum(isnan(data(Lside,:,i))) > 10
        data(Rside,:,i) = NaN;
        data(Lside,:,i) = NaN;
    else if abs(max(data(Rside,window,i))) < abs(max(data(Lside,window,i)))
            Rdata = data(Rside,:,i);
            data(Rside,:,i) = data(Lside,:,i);
            data(Lside,:,i) = Rdata;
            if iscell(alldata)
                Rsteplocs = steplocs(Rside,:,i);
                steplocs(Rside,:,i) = steplocs(Lside,:,i);
                steplocs(Lside,:,i) = Rsteplocs;
            end
        end
    end
end

sortedData = alldata;
if iscell(sortedData)
    sortedData{cond,1} = data;
    sortedData{cond,4} = steplocs;
else
    sortedData = data;
end
end


end