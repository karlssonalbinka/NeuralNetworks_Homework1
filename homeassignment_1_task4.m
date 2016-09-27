%homeassignment_1_task4
%Home assignment 1
%Task 4

%steps
%import training & validation data
% - normalize to zero mean and unit varianse (for both sets?!?)
%initialize weights randomly uniformly on [-0.2, 0.2]
%initialize thresholds randomly in [-1, 1]

% Make 100 independent training experiments,
% each with 2*10^5 iterations (one iteration corresponds to feeding a randomly
% chosen pattern and updating weights and thresholds). For each training ex-
% periment determine the minimum classification error for the training and the
% classification sets. Average these errors over the independent training exper-
% iments

% Train network with no hidden layers
% - assyncronous updating


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
% nbrIteration = 60000;
nbrIteration = 2*10^5;
% nbrExperiments = 2;
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
        
        w = w + lStep*delta_w;
        t = t - lStep*delta_t;
        
        %%%%%%%%%%%%%%%%%%%%
        %-------------- test Classification Error--------------------
%         classErr(nIteration) = sum(abs(tData(:,3) - sign(tanh(Beta*(tData(:,1:2)*w' - t))) ))/(2*length(tData));
        tmp = sum(abs(tData(:,3) - sign(tanh(Beta*(tData(:,1:2)*w' - t))) ))/(2*length(tData));
        if (tmp < minErr_t)
            minErr_t = tmp;
        end
        tmp = sum(abs(vData(:,3) - sign(tanh(Beta*(vData(:,1:2)*w' - t))) ))/(2*length(vData));
        if (tmp < minErr_v)
            minErr_v = tmp;
        end
        
        
%         tmp = 0;
%         for i = 1:length(tData)
%             b = w*tData(i,1:2)' - t;  % is this right with t (theta)?
%             o = tanh(Beta*b);
%             tmp = tmp + abs(tData(i,3)-sign(o));
%         end
%         classErr(nIteration) = tmp/(2*length(tData));
        %%%%%%%%%%%%%%%%%%%%
        
    end
    classErrMin_t(nExperiments) = minErr_t; %minimum classErr in training
    classErrMin_v(nExperiments) = minErr_v; %minimum classErr in validation
end
% plot(classErr)

% save('task4aResult', 'classErrMin_v', 'classErrMin_t');

%%

load task4aResult.mat

mean_t = mean(classErrMin_t)
mean_v = mean(classErrMin_v)
