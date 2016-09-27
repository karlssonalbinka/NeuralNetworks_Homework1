function error = calcClassificationError(patterns, outputs)
    %Takes patterns and their respective output and calculates the
    %classification error as mentioned in the task description for Example
    %sheet 1, task 4
    error = 0;
    nbrP = size(patterns,1);
    
    %take all zeros and make them positive 1
    outputs(outputs == 0) = 1;
    
    for i = 1:nbrP
        error = error + abs(patterns(i,3) - outputs(i,3)); %OBS may change depending on how outputs look
    end
    error = error/(2*nbrP);
end