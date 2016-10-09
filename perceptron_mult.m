% Perceptron for multiple classes

% input - X - matrix of examples in rows
%         Y - column vector of labels values are from 1 .. k
%         numLoops - max number of epochs (runs through all examples)

% output - Wout - weights for each class
%          D_X  - X data, with 1's appended to first column
%          D_Y  - Y data

function [Wout D_X D_Y] = perceptron_mult(X, Y, maxLoops)



X = [ ones(size(X,1),1), X];     % add in a column of ones for the bias term

numClasses = max(Y);             % picks max of Y vector, assumes max class number is present in Y



W = zeros(size(X,2),numClasses);  % weight matrix, number of columns equal to number of classes

stor_dotproducts = zeros(1,numClasses);     % storage variables for dotproducts, max Value and index
maxVal = 0;  
maxIdx = 0;


for j = 1:maxLoops               % j indicates number of loops elapsed
    
    changesCount = 0;            % count number of changes to Weights for an epoch
    
    for i=1:size(X,1)            % loop through examples
        stor_dotproducts = X(i,:)*W;
        
        [maxVal maxIdx] = max(stor_dotproducts);    % [0 0 0] max is first zero, Idx=1
        
        if (maxIdx ~= Y(i))
            W(:,maxIdx) = W(:,maxIdx) - X(i,:)';   % decrease the misclassified weight vector
            W(:,Y(i)) = W(:,Y(i)) + X(i,:)';       % increase the correct class weight vector
            
            changesCount = changesCount + 1;
        end
                       % if correct classification, do not update
    end % i loop, examples
    
    if (changesCount==0)   % if no weight changes after run through all data
        break;
    end

end % j loop, epochs

display(' ');display(['Training complete, number of epochs = ' num2str(j)]);

errors = 0; idx = zeros(1, size(X,1));
for i=1:size(X,1)                % count errors
    [m idx(i)] = max(X(i,:)*W);
    if idx(i) ~= Y(i)
        errors = errors+1;
    end
end               
display(['Errors : ' num2str(errors) ]); 


Wout = W;
D_X = X; D_Y = Y;          % store data for output, X has a one appended to front of every example



