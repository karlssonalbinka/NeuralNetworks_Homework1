%homeassignment_1_task4

clc
clear all

%Parameters (do not touch)
Beta = 0.5;
lStep = 0.01; %learning step

%import training and validation data
%each row is a pattern.
%Col 1 & 2 is input and col 3 is desired output
tData = load('train_data_2016.txt');
vData = load('valid_data_2016.txt');

% set mean of validation and training data to 0
var_col1 = std([tData(:,1); vData(:,1)]);
var_col2 = std([tData(:,1); vData(:,1)]);
col_1_shift = mean([tData(:,1); vData(:,1)]);
col_2_shift = mean([tData(:,2); vData(:,2)]);
%%%
tData(:,1) = (tData(:,1) - col_1_shift)/var_col1;
tData(:,2) = (tData(:,2) - col_2_shift)/var_col2;
vData(:,1) = (vData(:,1) - col_1_shift)/var_col1;
vData(:,2) = (vData(:,2) - col_2_shift)/var_col2;


%For loop parameters
nbrIteration = 2*10^5;
nbrExperiments = 100;

%test
classErr = zeros(1,nbrIteration);
classErrMin_t = zeros(1,nbrExperiments);
classErrMin_v = zeros(1,nbrExperiments);
test_t = zeros(1,nbrIteration);

for nExperiments = 1:nbrExperiments
    disp(nExperiments);
    %create random weights & thresholds
    w = rand(1,2)*0.4 - 0.2;
    t = rand(1,1)*2 - 1;
    minErr_t = 10^5;
    minErr_v = 10^5;
    
    for nIteration = 1:nbrIteration
        
        %Random what pattern to feed the system
        randPattern = floor(rand(1,1)* length(tData) + 1);
        
        b = w*tData(randPattern,1:2)' - t;
        
        Output = tanh(Beta*b);
        
        delta_t = Beta*(tData(randPattern,3) - Output)*(1-tanh(Beta*b)^2);
        delta_w = delta_t*tData(randPattern,1:2);
        
        %update wights and biases
        w = w + lStep*delta_w;
        t = t - lStep*delta_t;
        
        % check and save if new minimum classification error is found
        tmp = sum(abs(tData(:,3) - sign(tanh(Beta*(tData(:,1:2)*w' - t))) ))/(2*length(tData));
        if (tmp < minErr_t)
            minErr_t = tmp;
        end
        tmp = sum(abs(vData(:,3) - sign(tanh(Beta*(vData(:,1:2)*w' - t))) ))/(2*length(vData));
        if (tmp < minErr_v)
            minErr_v = tmp;
        end
    end
    classErrMin_t(nExperiments) = minErr_t; %minimum classErr in training
    classErrMin_v(nExperiments) = minErr_v; %minimum classErr in validation
end

save('task4aResult', 'classErrMin_v', 'classErrMin_t');

%%
clc
clear all
load task4aResult.mat

mean_t = mean(classErrMin_t)
mean_v = mean(classErrMin_v)
