      
%%%%%%%%%     %%%   %%%%%     %%%   %%%%%     %%%   %%%   %%%%%     %%%     %%%%%%%%%%
%%%    %%%    %%%   %%%%%%    %%%   %%%%%%    %%%   %%%   %%%%%%    %%%    %%%       %%
%%%   %%%     %%%   %%% %%%   %%%   %%% %%%   %%%   %%%   %%% %%%   %%%   %%%
%%%  %%%      %%%   %%%  %%%  %%%   %%%  %%%  %%%   %%%   %%%  %%%  %%%   %%%
%%%    %%%    %%%   %%%   %%% %%%   %%%   %%% %%%   %%%   %%%   %%% %%%   %%%     %%%%%%
%%%     %%%   %%%   %%%    %%%%%%   %%%    %%%%%%   %%%   %%%    %%%%%%   %%%        %%% 
%%%    %%%    %%%   %%%     %%%%%   %%%     %%%%%   %%%   %%%     %%%%%    %%%      %%%
%%%%%%%%%     %%%   %%%      %%%%   %%%      %%%%   %%%   %%%      %%%%     %%%%%%%%%%


% with binary sort


%%% --- INPUT : row vector
%%% --- OUTPUT : row vector with values binned into 10 bins

function [x_2] = binning_12(input_vec)


% examples
% x = [-5, -4.4, -3.2, -1.2, -0.6, 0.1, 0.3, 0.65, 0.8, 1, 1.2, 3.3];
% 
% x = rand(1,10);
% x = zscore(x);

% plot(x, [0 0 0 0 0 0 0 0 0 0 0 0],'ko');


% 12 intervals [-inf, -1.5] [-1.5, -1.2] [-1.2, -0.9] [-0.9, -0.6] [-0.6, -0.3] [-0.3, 0]
%                [0, 0.3] [0.3, 0.6] [0.6, 0.9] [0.9, 1.2] [1.2, 1.5] [1.5, inf]

% 11 'boundaries'


% test a value ,       log(10), ~3 evaluations

% val = -inf

% examples in row vector
x = input_vec; 
if (size(x,1) > size(x,2)) % test for column vector
    x = x';
end

x_2=zeros(1,size(x,2));

for i = 1:size(x,2)
    val = x(i);
    

if (val >= 0)
    % display('larger than 0');
    if (val >= 0.9)
        % display('larger than .9');
        if (val >= 1.2)
            % display('larger than 1.2');
            if (val >= 1.5)
                % display('larger than 1.5');
                x_2(i)=1.5;
            else
                % display('smaller than 1.5');
                x_2(i)=1.2;
            end
        else
            % display('smaller than 1.2');
            x_2(i)=0.9;
        end
    else
        % display('smaller than 0.9');
        if (val >= 0.3)
            % display('larger than 0.3');
            if (val >= 0.6)
                % display('larger than 0.6');
                x_2(i)=0.6;
            else
                % display('smaller than 0.6');
                x_2(i)=0.3;
            end
                
        else
            % display('smaller than 0.3');
            x_2(i)=0;
        end
    end
        
else
    % display('smaller than 0');
    if (val >= -0.9)
        % display('larger than -0.9');
        if (val >= -0.6)
            % display('larger than -0.6');
            if (val >= -0.3)
                % display('larger than -0.3');
                x_2(i)=-0.3;
            else
                % display('smaller than -0.3');
                x_2(i)=-0.6;
            end
                
        else
            % display('smaller than -0.6');
            x_2(i)=-0.9;
        end
    else
        % display('smaller than -0.9');
        if (val >= -1.2)
            % display('larger than -1.2');
            x_2(i)=-1.2;
        else
            % display('smaller than -1.2');
            x_2(i)=-1.5;
        end
    end
end


end %for loop, going through x's elements 




