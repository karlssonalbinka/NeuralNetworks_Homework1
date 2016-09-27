%homeassignment_1_task4
%Home assignment 1
%Task 4 b

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
neurons = [2 4 8 16 32];

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
nbrIteration = 30000;
% nbrIteration = 2*10^5;
nbrExperiments = 2;
nbrNeurons = 1;
% nbrExperiments = 100;

%test
classErrMin_t = zeros(5,nbrExperiments);
classErrMin_v = zeros(5,nbrExperiments);

for nNeurons = 1:nbrNeurons
    
    for nExperiments = 1:nbrExperiments
        disp(nExperiments);
        %create random weights & thresholds
        w1 = rand(2,neurons(nNeurons) )*0.4 - 0.2;
        w2 = rand(1,neurons(nNeurons) )*0.4 - 0.2;
        t1 = rand(2,1)*2 - 1;
        t2 = rand(1,1)*2 - 1;
        minErr_t = 10^5;
        minErr_v = 10^5;
        
        for nIteration = 1:nbrIteration
            
            %Random what pattern to feed the system
            randPattern = floor(rand(1,1)* length(tData) + 1)
            disp('1')
            b1 = w1*tData(randPattern,1:2)' - t1
            disp('2')
            V = tanh(Beta*b1) %The output to the hidden layer
            disp('3')
            b2 = w2*V - t2
            disp('4')
            Output = tanh(Beta*b2) %The output to the hidden layer
            
            
            disp('5')
            delta_t2 = Beta*(tData(randPattern,3) - Output)*(1-tanh(Beta*b2)^2);
            disp('6')
            delta_w2 = delta_t2*V;
            disp('7')
            delta_t1 = Beta*w2*delta_w2*(1-tanh(Beta*b1).^2);
            disp('8')
            delta_w1 = delta_t1*tData(randPattern,1:2);
            
            disp('9')
            w1 = w1 + lStep*delta_w1;
            disp('10')
            t1 = t1 - lStep*delta_t1;
            disp('11')
            w2 = w2 + lStep*delta_w2';
            disp('12')
            t2 = t2 - lStep*delta_t2;
            
            %%%%%%%%%%%%%%%%%%%%
            %-------------- test Classification Error--------------------
            b1 = w1*tData(:,1:2)';
            disp('20')
            b1 = [b1(1,:)-t1(1); b1(2,:)-t1(2)];
            disp('21')
            V = tanh(Beta*b1); %The output to the hidden layer
            disp('22')
            b2 = w2*V - t2;
            disp('23')
            Output = tanh(Beta*b2); %The output to the hidden layer
            disp('24')
            % CONTROLLERA 
%             classErr(nIteration) = sum(abs(tData(:,3) - sign(tanh(Beta*(tData(:,1:2)*w' - t))) ))/(2*length(tData));
            tmp = sum(abs(tData(:,3) - sign( tanh(Beta*b2) )))/(2*length(tData));
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
end
% plot(classErr)

save('task4aResult', 'classErrMin_v', 'classErrMin_t');

%%
% clear all
% clc
% neurons = [2 4 8 16 32];
% for i = 1:2
%     nbr = 2+neurons(i)+1;
%     w = rand(nbr,nbr )*0.4 - 0.2;
%     w = triu(w,1);
%     w(1:2,
% end
