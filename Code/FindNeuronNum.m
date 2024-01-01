function neuron_num = FindNeuronNum(currneur, relevant_neurons)
neurons{1} = currneur;
num_of_sides = 1;
if length(currneur) == 3
    neurons{1} = strcat(currneur,'R');
    neurons{2} = strcat(currneur,'L');
    num_of_sides = 2;
end

for i = 1:num_of_sides
indx = strfind(relevant_neurons, neurons{i});
neuron_num(i) = find(not(cellfun('isempty',indx)));

end


 

end