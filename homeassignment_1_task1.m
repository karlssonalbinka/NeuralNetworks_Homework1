% task 1
clc
clear all;
close all

%changeable parameters
nbrTimes = 1000;  %nbr experiments to average over

%parameters
N = [100, 200]; %nbr of bits
P = [10 20 30 40 50 75 100 150 200]; %nbr of patterns
loopLength = ceil(nbrTimes./P); %To ensure that all patterns get an average of 1000 experiments

percentWrong = zeros(length(N), length(P) );    
alpha = zeros(length(N), length(P));

for N_ITERATION = 1:length(N)
    for P_ITERATION = 1:length(P)
        
        totBits = N(N_ITERATION)*P(P_ITERATION);
        alpha(N_ITERATION, P_ITERATION) = P( P_ITERATION )/N( N_ITERATION);
        tot = 0;
        for ITERATION = 1:loopLength(P_ITERATION)
            
            %Create random patterns:
            patterns = sign( round( rand(P( P_ITERATION ), N( N_ITERATION) )) - 0.1 );
            
            
            %Calculate weights
            w = patterns'*patterns/N(N_ITERATION);
            for i = 1: length(w)
                w(i,i) = 0;
            end
            nextPatterns = patterns*w;
            nextPatterns( nextPatterns == 0) = 1; %take care of elements that equal 0
            nextPatterns = sign( nextPatterns );
            
            wrongBits = totBits - sum( sum(nextPatterns == patterns),2);
            tot = tot+wrongBits;
        end
            percentWrong(N_ITERATION, P_ITERATION) = tot/(totBits*loopLength(P_ITERATION)); %get mean of wrong percent
    end
end

%Create throretical P_error
X = [0.01:0.02:3];
Y= zeros(1, length(X));
for i = 1:length(X);
    Y(i) = 1-erf( sqrt(1/(2*X(i))) ); %1/X = N/P
end
Y = Y/2;
%%
hold on
plot(X, Y*100, 'k');
plot(alpha(1,:), percentWrong(1,:)*100, '-*r');
plot(alpha(2,:), percentWrong(2,:)*100, '-o');
title('One step error estimate using the Hopfield model');
ylabel('Percent Error (%)');
xlabel('\alpha (p/N)');
legend('Theoretical P_{error}', 'N=100', 'N=200', 'Location', 'northwest');


