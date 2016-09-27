clc
clear all;
colormap gray

% zero = imread('zero.jpg');
% one = imread('one.jpg');
% two = imread('two.jpg');
% three = imread('three.jpg');
% four = imread('four.jpg');
zero = imread('0.gif');
one = imread('1.gif');
two = imread('2.gif');
three = imread('3.gif');
four = imread('4.gif');

%making it a grayscale (2D) image
zero = im2bw(zero);
one = im2bw(one);
two = im2bw(two);
three = im2bw(three);
four = im2bw(four);

tmp = {zero, one, two, three, four};

%parameters
nbrRepeat = 100;
nbrP = 5;       % number of patterns
maxGen = 5;    % the maximum number of generations
switchedBits = 120;  %maximum number of switched bits
[rows, cols] = size(zero);  %dimension of our figures
N = rows*cols;  %nbr of bits
percent = zeros(nbrP,switchedBits); %the stored percentages

%remake to vektors:
pattern = zeros(nbrP, N);
for i=1:nbrP
    for j=1:rows
        for k=1:cols
            pattern(i, (j-1)*cols + k) = tmp{i}(j, k);
        end
    end
end
%convert to +/-1 in stead of 1,0
pattern = sign(pattern-0.1);

%Calculate weights
w = pattern'*pattern/N;


for PATTERN_NBR = 1:nbrP
    PATTERN_NBR
    for nbrSwitchedBits=1:switchedBits
        %         disp(strcat('   - ',  num2str(nbrSwitchedBits)) );
        
        for k = 1:nbrRepeat
            %get distorted patterns
            switchBit = floor(rand(nbrRepeat,nbrSwitchedBits)*N) + 1;  %OBS it might switch the same twice. Fix later
            for i=1:nbrRepeat
                tmp = pattern(PATTERN_NBR,:);
                tmp(switchBit(i,:)) = tmp( switchBit(i,:) )*(-1);
                distP(i, :) = tmp;
                %diffPix = N - sum(sum(distP(i,:) == pattern(PATTERN_NBR,:)),2)
            end
            
            % Draw distorted pattern before change
%             figure = zeros(rows, cols);
%             for i=1:rows
%                 index = (i-1)*cols;
%                 figure(i,:) = distP( k, index+1:index+cols);
%             end
%             figure = figure == 1;   %remake the figure to 0, 1 values
%             image(figure, 'CDataMapping','scaled');
%             drawnow
%             pause(1);
            
            for j = 1:maxGen;
                sequence = randperm(N);
                for i = 1:N
                    tmp = w(sequence(i),:)*distP(k,:)';
                    if ( tmp ~= 0 )
                        distP(k,sequence(i)) = sign(tmp);
                    else
                        distP(k,sequence(i)) = 1;
                    end
                end
                
                
                % Remake figure for understanding
%                 figure = zeros(rows, cols);
%                 for i=1:rows
%                     index = (i-1)*cols;
%                     figure(i,:) = distP( k, index+1:index+cols);
%                 end
%                 figure = figure == 1;   %remake the figure to 0, 1 values
%                 image(figure, 'CDataMapping','scaled');
%                 drawnow
%                 pause(.15);
            end
            missmatchedPixels = N - sum(sum(distP(k,:) == pattern(PATTERN_NBR,:)),2);
            if ( missmatchedPixels == 0 )
                percent(PATTERN_NBR, nbrSwitchedBits) = percent(PATTERN_NBR, nbrSwitchedBits)+ 1;
            end
        end
        
    end
end

percent = percent/nbrRepeat;
plot(percent');
legend('Zero', 'One', 'Two', 'Three', 'Four');
xlabel('Number flipped bits');
ylabel('Percent of correct matches (%)');
%did_it_work = sum(sum(distP == pattern(1,:)),2)
