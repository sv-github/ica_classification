% Training ICA multi class classifier

%     --- Naive Bayes ---  with 15x15x15 image matrix resizing

%%% NOTE: requires viewer labeled training data, components picked out from
%%% GIFT_ICA output that are auditory, visual, default-mode, motor


% load 'libraries' (nifti_lab)
% addpath(genpath('I:\VPGroup\Data\fMRI_Studies\StrokePlasticity\Scans\svyat\NIFTI_mlab'));
addpath(genpath('./NIFTI_mlab'));


% NOTE: data loading is commented out. the naive bayes algorithm starts on
% line 180 below

%{
% set training subjects (10 here as an example) and testing subjects

% set up data matrix X and label matrix Y for training          
% 191 components, 3375 (15*15*15 attributes)
X = zeros(191, 15*15*15);   % <----- change   191 to number of components in training set
Y = zeros(191,1);           % <----- change   191 to number of components in training set



% --- load training examples

% --- load Visual components
%  S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S
% S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S


sub_nums = [1 2 3 4 5 6 7 8 9 10];   % training subject list
                                            
%cd('to the directory with manually set up component matrix ');  
load('visual_ic_matrix.mat');                % manually filled in
                                             % variable = ic_matrix_vis,  components in columns for each subject
                                             
%                                          subj1 [IC1 IC2 IC15 IC16 IC27]      [1 2 15 16 27]
%                                          subj2 [IC2 IC4 IC15 0    0   ]      [2 4 15  0  0]
%                                          subj3 [IC4 IC6 IC7  0    0   ]  --> [4 6  7  0  0]
%                                          .                                    .
%                                          .                                    .
%                                          .     [                      ]      [5 8 11 20  0]
%                                                these are all components
%                                                that are visual network

%                                        NOTE: an example visual_ic_matrix
%                                        is provided
comp_count_vis = 0;

for i = 1:size(sub_nums,2)             % example directory where each subject has 
                                       % the GIFT ICA output nifti file of components
                                       % called 'gift_sub01_component_ica_s1_.ni
                                       
    %sub_dir = sprintf('C:\\Users\\vxs207\\Documents\\GIFT_23sub\\sub%g',sub_nums(i));
    cd(sub_dir);   % change to subject directory
    load_nii('gift_sub01_component_ica_s1_.nii');  ic = ans.img; ic = double(ic); % load components
    
    for j = 1:size(find(ic_matrix_vis(i,:)),2)   % count nonzero elements only
            comp_count_vis = comp_count_vis+1;
            
            % load component, get 15x15x15 input
            comp = ic(:,:,:,ic_matrix_vis(i,j));
            comp_rs = imresize3d(comp, 15); comp_rs_lin = comp_rs(:);   % resize component maps
            X(comp_count_vis,:) = comp_rs_lin;            Y(comp_count_vis) = 1;     %visual
            
        
    end
    
end      % comp_count_vis = 43,     43 so far (for example)



% --- load Sensorimotor components
%  S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S
% S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S


                                          
%cd('to the directory with manually set up component matrix ');  
load('motor_ic_matrix.mat');  % manually filled in as before
        % variable = ic_matrix                   components in columns for each subject
        
comp_count_motor = 0;

for i = 1:size(sub_nums,2)
    %sub_dir = sprintf('C:\\Users\\vxs207\\Documents\\GIFT_23sub\\sub%g',sub_nums(i));
    cd(sub_dir);
    load_nii('gift_sub01_component_ica_s1_.nii');  ic = ans.img; ic = double(ic); % load components
    
    for j = 1:size(find(ic_matrix(i,:)),2)   % count nonzero elements only
            comp_count_motor = comp_count_motor+1;
            
            % load component
            comp = ic(:,:,:,ic_matrix(i,j));
            comp_rs = imresize3d(comp, 15); comp_rs_lin = comp_rs(:);
            X(comp_count_vis+comp_count_motor,:) = comp_rs_lin;            Y(comp_count_vis+comp_count_motor) = 2;     %motor
            
        
    end
    
end      % comp_count = 30,     43+30 = 73 components so far (for example)

                                         
                                  

% --- load Default mode components
%  S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S
% S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S
                                           


%cd('to the directory with manually set up component matrix ');   
load('defmode_ic_matrix.mat');  % manually filled in as before
   % variable = ic_matrix_defmode

comp_count_defmode = 0;

for i = 1:size(sub_nums,2)
    %sub_dir = sprintf('C:\\Users\\vxs207\\Documents\\GIFT_23sub\\sub%g',sub_nums(i));
    cd(sub_dir);
    load_nii('gift_sub01_component_ica_s1_.nii');  ic = ans.img; ic = double(ic); % load components
    
    for j = 1:size(find(ic_matrix_defmode(i,:)),2)   % count nonzero elements only
            comp_count_defmode = comp_count_defmode+1;
            
            % load component
            comp = ic(:,:,:,ic_matrix_defmode(i,j));
            comp_rs = imresize3d(comp, 15); comp_rs_lin = comp_rs(:);
            X(comp_count_vis+comp_count_motor+comp_count_defmode,:) = comp_rs_lin;   Y(comp_count_vis+comp_count_motor+comp_count_defmode) = 3; % default mode
            
        
    end
    
end      % comp_count_defmode = 85,     43+30+85 = 158 components (for example)




% --- load Auditory components
%  S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S
% S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S S



%cd('to the directory with manually set up component matrix ');
load('auditory_ic_matrix.mat');  % manually filled in as before
     % variable = ic_matrix_auditory

comp_count_auditory = 0;

for i = 1:size(sub_nums,2)
    %sub_dir = sprintf('C:\\Users\\vxs207\\Documents\\GIFT_23sub\\sub%g',sub_nums(i));
    cd(sub_dir);
    load_nii('gift_sub01_component_ica_s1_.nii');  ic = ans.img; ic = double(ic); % load components
    
    for j = 1:size(find(ic_matrix_auditory(i,:)),2)   % count nonzero elements only
            comp_count_auditory = comp_count_auditory+1;
            
            % load component
            comp = ic(:,:,:,ic_matrix_auditory(i,j));
            comp_rs = imresize3d(comp, 15); comp_rs_lin = comp_rs(:);
            X(comp_count_vis+comp_count_motor+comp_count_defmode+comp_count_auditory,:) = comp_rs_lin;   Y(comp_count_vis+comp_count_motor+comp_count_defmode+comp_count_auditory) = 4; % auditory
            
        
    end
    
end      % comp_count_auditory = 33,     43+30+85+33 = 191 components (for example)


% --- Save training data ---
cd('data directory');          % change to data directory (can be where subjects are saved)
save('IC_training_data_X.mat', 'X');
save('IC_training_data_Y.mat', 'Y');
%}


%%% --- prepare data for naive bayes algorithm (requires data to be in discrete intervals)
%%% --- standardize nonzero elements of examples



%%% --- bin data using specified intervals (with binning_12.m helper function)
%       convert value to a bin value

%X_1b=binning_12(X_1);

%x_xb = [X_1; X_1b];  matrix can be viewed for comparison

% --- convert all examples to binned discrete values
X_1b = zeros(size(X,1),size(X,2));
X_1 = zeros(1,size(X,2));

for i = 1:size(X,1)
    Nz = X(i,:)~=0;      % only look at nonzero values (look only inside the brain)
    X_1 = X(i,:);       % copy example over
    X_1(Nz) = zscore(X_1(Nz));  % zscore only the nonzero values
    
    X_1b(i,:) = binning_12(X_1);  % X_1binned
    
    clear Nz;
end




%%% --- 'convert' attribute values to integers (1..11)
X_1b = (((X_1b + 1.5).*(10/3))+1);        % (-1.5 + 1.5)*(10/3) -> 0 
                                          % (-1.2 + 1.5)*(10/3) -> 1 

test_int = uint8(X_1b(1,1));

X_1b = uint8(X_1b);

% initialize counter variables
class1_total= 0;       
class2_total= 0; class3_total= 0; class4_total= 0;    % count of number of examples in each class

attrib_counter_c1 = zeros(11,size(X_1b,2)); attrib_counter_c2 = zeros(11,size(X_1b,2)); 
attrib_counter_c3 = zeros(11,size(X_1b,2)); attrib_counter_c4 = zeros(11,size(X_1b,2));
val=0;

%%% --- loop through components (examples) and attributes, count up frequencies of values
for cmpts = 1:size(X_1b,1) 
    
    if Y(cmpts)==1        % class 1 - visual
        class1_total = class1_total+1;
    
        for att = 1:size(X_1b,2)
            val = X_1b(cmpts,att);
            attrib_counter_c1(val,att) = attrib_counter_c1(val,att)+1;
        end
    
    elseif Y(cmpts)==2    % class 2 - motor
        class2_total = class2_total+1;
        
        for att = 1:size(X_1b,2)
            val = X_1b(cmpts,att);
            attrib_counter_c2(val,att) = attrib_counter_c2(val,att)+1;
        end
             
    elseif Y(cmpts)==3    % class 3 - default mode
        class3_total = class3_total+1;
        
        for att = 1:size(X_1b,2)
            val = X_1b(cmpts,att);
            attrib_counter_c3(val,att) = attrib_counter_c3(val,att)+1;
        end
        
    else                  % class 4 - auditory
        class4_total = class4_total+1;
        
        for att = 1:size(X_1b,2)
            val = X_1b(cmpts,att);
            attrib_counter_c4(val,att) = attrib_counter_c4(val,att)+1;
        end
        
    end %if,else switch for 4 classes 
    
end %forloop
        


%%% --- Log probability summing     (class_out = argmax cj in C  P(cj) PRODCUT(P(ai|cj)))

% compute probability tables from attrib_counter tables
attrib_prob_c1=zeros(11,size(X_1b,2)); attrib_prob_c2=zeros(11,size(X_1b,2));
attrib_prob_c3=zeros(11,size(X_1b,2)); attrib_prob_c4=zeros(11,size(X_1b,2));

attrib_prob_c1 = (attrib_counter_c1 + 0.25)/(class1_total + 1);
attrib_prob_c2 = (attrib_counter_c2 + 0.25)/(class2_total + 1);
attrib_prob_c3 = (attrib_counter_c3 + 0.25)/(class3_total + 1);
attrib_prob_c4 = (attrib_counter_c4 + 0.25)/(class4_total + 1);

Prob_c1 = class1_total/cmpts; Prob_c2 = class2_total/cmpts;
Prob_c3 = class3_total/cmpts; Prob_c4 = class4_total/cmpts;






% ____Testing_____  (testing IC maps, using imresize3d.m and sum_logs_4class.m helper functions)
%                   requires processed data and ICA maps

%cd('testing subject directory');
    
load_nii('gift_sub01_component_ica_s1_.nii');  ic = ans.img; ic = double(ic); % load components


% --- Load testing subject

comp_count_vis = 0;
                            
vis_components = [1 5 17];  % example visual component numbers
                 % NOTE: update to subject's specific component numbers
                 % OR can create a component matrix as for training subjects

for j = 1:size(vis_components,2)   
        comp_count_vis = comp_count_vis+1;

        % load component, get 15x15x15 input
        comp = ic(:,:,:,vis_components(j));
        comp_rs = imresize3d(comp, 15); comp_rs_lin = comp_rs(:);   % resize component maps
        X_leftOut(comp_count_vis,:) = comp_rs_lin;            Y_leftOut(comp_count_vis) = 1;     %visual


end


comp_count_motor = 0;
                            
motor_components = [3 15 27];  % example visual component numbers
                 % NOTE: update to subject's specific component numbers

for j = 1:size(motor_components,2)   
        comp_count_motor = comp_count_motor+1;

        % load component, get 15x15x15 input
        comp = ic(:,:,:,motor_components(j));
        comp_rs = imresize3d(comp, 15); comp_rs_lin = comp_rs(:);   % resize component maps
        X_leftOut(comp_count_vis+comp_count_motor,:) = comp_rs_lin;            Y_leftOut(comp_count_vis+comp_count_motor) = 2;     %motor


end


comp_count_defmode = 0;
                            
defmode_components = [9 12 20];  % example visual component numbers
                 % NOTE: update to subject's specific component numbers

for j = 1:size(defmode_components,2)   
        comp_count_defmode = comp_count_defmode+1;

        % load component, get 15x15x15 input
        comp = ic(:,:,:,defmode_components(j));
        comp_rs = imresize3d(comp, 15); comp_rs_lin = comp_rs(:);   % resize component maps
        X_leftOut(comp_count_vis+comp_count_motor+comp_count_defmode,:) = comp_rs_lin;            Y_leftOut(comp_count_vis+comp_count_motor+comp_count_defmode) = 3;     %default-mode


end


comp_count_auditory = 0;
                            
auditory_components = [10 18 19];  % example visual component numbers
                 % NOTE: update to subject's specific component numbers

for j = 1:size(auditory_components,2)   
        comp_count_auditory = comp_count_auditory+1;

        % load component, get 15x15x15 input
        comp = ic(:,:,:,auditory_components(j));
        comp_rs = imresize3d(comp, 15); comp_rs_lin = comp_rs(:);   % resize component maps
        X_leftOut(comp_count_vis+comp_count_motor+comp_count_defmode+comp_count_auditory,:) = comp_rs_lin;            Y_leftOut(comp_count_vis+comp_count_motor+comp_count_defmode+comp_count_auditory) = 4;     %auditory


end


% --- save testing subject data ---
save('IC_testing_data_X.mat', 'X_leftOut');
save('IC_testing_data_Y.mat', 'Y_leftOut');




% variable X_leftOut  15*15*15 elements in each example
Y_leftOut=Y_leftOut'; % variable Y_leftOut, column vector  in 1,2,3,4 form

% --- rename X,Y_leftOut to match variables below
testing_X = X_leftOut;    testing_Y = Y_leftOut;    clear X_leftOut;  clear Y_leftOut;

label_output = zeros(size(testing_Y));   % variable to hold classification output



% --- zscore and binning on testing_X
% --- convert all examples to binned discrete values
testing_X_1b = zeros(size(testing_X,1),size(testing_X,2));
testing_X_1 = zeros(1,size(testing_X,2));

for i = 1:size(testing_X,1)
    t_Nz = testing_X(i,:)~=0;      % only look at nonzero values (look only inside the brain)
    testing_X_1 = testing_X(i,:);       % copy example over
    testing_X_1(t_Nz) = zscore(testing_X_1(t_Nz));  % zscore only the nonzero values

    testing_X_1b(i,:) = binning_12(testing_X_1);  % testing_X_1binned

    clear t_Nz;
end


%%% --- 'convert' attribute values to integers (1..11)
testing_X_1b = (((testing_X_1b + 1.5).*(10/3))+1);        % (-1.5 + 1.5)*(10/3) -> 0
% (-1.2 + 1.5)*(10/3) -> 1

testing_X_1b = uint8(testing_X_1b);




%%% --- test the 4 classes(networks)
% using label_output
for i = 1:size(testing_X_1b,1)
    [max_class_o idx_class_o] = sum_logs_4class(testing_X_1b(i,:),Prob_c1,Prob_c2,Prob_c3,Prob_c4,attrib_prob_c1,attrib_prob_c2,attrib_prob_c3,attrib_prob_c4);
    label_output(i) = idx_class_o;
    clear max_class_o; clear idx_class_o;    
end   


error= (sum(testing_Y ~= label_output)/size(testing_Y,1));

display('error of classification'); % human labels - testing_Y,   algorithm labels - label_output 
display(error);

