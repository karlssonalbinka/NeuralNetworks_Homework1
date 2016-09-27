% task 3a
clc
clear all;
close all


%changeable parameters
nbrTimes = 100;
nbrGenerations = 200;
N = 500;
B = 2;
P = [10 20 30 40 50 75 100 150 200 300 400 500];
% P=20; %test
order = zeros(1, length(P));
g = @(b) 1/(exp(-2*B*b) + 1);   % Probability function

%parameters
alpha = zeros(length(N), length(P));

for N_ITERATION = 1:length(N) %will be 1:1 in task 3 a)
    nbrBits = N( N_ITERATION);
    for P_ITERATION = 1:length(P)
        nbrPatterns = P( P_ITERATION );
        
%         totBits = nbrBits*nbrPatterns;
        alpha(N_ITERATION, P_ITERATION) = nbrPatterns/nbrBits;

        nbrRepeat = round(nbrTimes/nbrPatterns); %To get a mean out of nbrTimes times (because the number of Patterns change that we take mean of
        if nbrRepeat < 1    % if nbrPatters is > 2*nbrTimes
            nbrRepeat = 1;
        end
        for REPEAT_ITERATION = 1:nbrRepeat
            
            %Create random patterns:
            patterns = sign( round( rand(P( P_ITERATION ), N( N_ITERATION) )) - 0.1 );
            
            %Calculate weights
            w = patterns'*patterns/nbrBits;
            for i = 1: length(w)
                w(i,i) = 0;
            end
            
            %Perform iterations
            updatedPattern = patterns;
            for i = 1:nbrGenerations
                
                %get the order in which to update the bits:
                sequence = randperm(N);
                randNr = rand(nbrPatterns,N);
                for j = 1:N
                    tmp = updatedPattern*w(sequence(j) ,:)';
                    
                    %Stochastic part
                    for k = 1:nbrPatterns
                        if ( g(tmp(k)) > randNr(k,j) )
                            updatedPattern(k, sequence(j)) = 1;
                        else
                            updatedPattern(k, sequence(j)) = -1;
                        end
                    end
                end
                
                if(i > 150) %make a better statement later
                    %Store order parameter for each generation
                    tmp = 0;
                    for j = 1:nbrPatterns
                        tmp = tmp + updatedPattern(j,:)*patterns(j,:)';
                    end
                    m(i-150) = tmp/(nbrPatterns*nbrBits);
                end
%                 m(i) = tmp/(nbrPatterns*nbrBits);
                
            end %nbrGenerations
            order(P_ITERATION) = order(P_ITERATION) + mean(m);
        end %REPEAT_ITERATION
        order(P_ITERATION) = order(P_ITERATION)/nbrRepeat;
    end %P_ITERATION
end %N_ITERATION

%save('orderAndAlpha_task3a.mat', 'order', 'alpha');

hold on
plot(alpha, order);
legend('order parameter');
xlabel('\alpha');
ylabel('order parameter');

%% Load and plot the saved result from previous section
load orderAndAlpha_task3a.mat


plot(alpha, order);
legend('m(\beta=2)', 'Interpreter', 'latex'); %doesn't work.
% set(h, 'Interpreter', 'latex')
xlabel('\alpha');
ylabel('order parameter');

%%  3b
clc
clear all;
close all


%changeable parameters
nbrTimes = 100;
nbrGenerations = 200;
N = [50 100 250];
B = 2;

%Get same alpha as before
oldP = [10 20 30 40 50 75 100 150 200 300 400 500];
percent=oldP/500;
P = [percent*N(1); percent*N(2); percent*N(3)];
P = round(P);

%parameters
% P=20; %test

alpha = zeros(length(N), length(P));
order = zeros(length(N), length(P));
g = @(b) 1/(exp(-2*B*b) + 1);   % Probability function


for N_ITERATION = 1:length(N)
    disp(['N = ' num2str(N(N_ITERATION) )])
    nbrBits = N( N_ITERATION);
    for P_ITERATION = 1:length(P)
        nbrPatterns = P(N_ITERATION, P_ITERATION );
        
        alpha(N_ITERATION, P_ITERATION) = nbrPatterns/nbrBits;

        nbrRepeat = round(nbrTimes/nbrPatterns); %To get a mean out of nbrTimes times (because the number of Patterns change that we take mean of
        if nbrRepeat < 1    % if nbrPatters is > 2*nbrTimes
            nbrRepeat = 1;
        end
        for REPEAT_ITERATION = 1: nbrRepeat 
            
            %Create random patterns:
            patterns = sign( round( rand(P(N_ITERATION, P_ITERATION ), N( N_ITERATION) )) - 0.1 );

            %Calculate weights
            w = patterns'*patterns/nbrBits;
            for i = 1: length(w)
                w(i,i) = 0;
            end

            %Perform iterations
            updatedPattern = patterns;
            for i = 1:nbrGenerations

                %get the order in which to update the bits:
                sequence = randperm(N(N_ITERATION));
                randNr = rand(nbrPatterns,N(N_ITERATION));
                for j = 1:N(N_ITERATION)
                    tmp = updatedPattern*w(sequence(j) ,:)';
                    %Stochastic part
                    for k = 1:nbrPatterns
                        if ( g(tmp(k)) > randNr(k,j) )
                            updatedPattern(k, sequence(j)) = 1;
                        else
                            updatedPattern(k, sequence(j)) = -1;
                        end
                    end
                end
                
                if(i > 150) %make a better statement later
                    %Store order parameter for each generation
                    tmp = 0;
                    for j = 1:nbrPatterns
                        tmp = tmp + updatedPattern(j,:)*patterns(j,:)';
                    end
                    m(i-150) = tmp/(nbrPatterns*nbrBits);
                end
%                 m(i) = tmp/(nbrPatterns*nbrBits);
                
            end %nbrGenerations
            
            order(N_ITERATION, P_ITERATION) = order(N_ITERATION, P_ITERATION) + mean(m);
        end %REPEAT_ITERATION
        order(N_ITERATION, P_ITERATION) = order(N_ITERATION, P_ITERATION)/nbrRepeat;
    end %P_ITERATION
end %N_ITERATION

save('orderAndAlpha_task3b.mat', 'alpha', 'order');

hold on
plot(alpha', order');
legend('N = 50', 'N = 100', 'N = 250');
xlabel('\alpha');
ylabel('order parameter');

%% Load and plot result from 3b
load orderAndAlpha_task3b.mat;
plot(alpha', order');
legend('N = 50', 'N = 100', 'N = 250');
xlabel('\alpha');
ylabel('order parameter');

%% Load and plot everything in same
%% Load and plot the saved result from previous section
clc
clear all
t1 = load('orderAndAlpha_task3a.mat');
t2 = load ('orderAndAlpha_task3b.mat');
order = [t2.order; t1.order];
alpha = [t2.alpha; t1.alpha];

hold on
plot(alpha(1,:)', order(1,:)', '-');
plot(alpha(2,:)', order(2,:)', ':g');
plot(alpha(3,:)', order(3,:)', '*k');
plot(alpha(4,:)', order(4,:)', '-*r');

legend('N = 50', 'N = 100', 'N = 250', 'N = 500'); %doesn't work.
% set(h, 'Interpreter', 'latex')
xlabel('\alpha');
ylabel('order parameter');