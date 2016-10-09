                  %%%%% Processing %%%%%

% Creates blurred, 3dresampled resting data to be used in GIFT ICA
% extraction

% *required resting data files:      rest.nii  ,    anat.nii
% **requires MATLAB to be run on linux computer with AFNI and FSL installed
    
% INPUT: rest.nii, anat.nii
% OUTPUT: master_cut.nii, GIFT_ICA maps


% --- rest and anat files needed for processing below VVVVVVVVV

          % --- change to subject directory
%cd(subject_directory);

mkdir('proc_Calh_templates');       % make processing directory
          

if ~exist('anat+orig.BRIK', 'file')     % create anat+orig if does not exist
    !3dcopy anat.nii anat
end
    

% --- copy resting file and anat file into proc_Calh_templates
!mv rest.nii ./proc_Calh_templates
%!mv rest+orig* ./proc_Calh_templates
!mv anat+orig* ./proc_Calh_templates
                 % anat has skull
                 
cd('./proc_Calh_templates');        % switch to processing dir, rename files
                       

             % --- begin processing of rest file and templates ---
             
%%% --- NOTE: RPI rotate before beginning processing --- %%%

% ( 1 ) cut first 4 volumes----------------------cut first 4 volumes ---------------
!3dcalc -a rest.nii[4..$] -expr 'a' -prefix rest_x.nii
                            
    
% ( 2 ) motion corr--------------------------- & slice time correction (-tshift)------
!3dvolreg -prefix rest_tshift.nii -tshift -Fourier -verbose -base rest_x.nii[0] -zpad 4 -dfile rest_motion.txt rest_x.nii
% !1dplot -volreg rest_motion.txt[1..6]
                           %3dhistog - histogram

                           
% ( 3 ) skull strip rest file ---------------------------------------------------
!bet rest_tshift.nii rest_ss.nii -F -f 0.5 -g 0


% prepare file for transformation, normalization
!gunzip rest_ss.nii.gz
!3dcopy rest_ss.nii rest_ss

                       

% ( 4 ) normalize  anat to standard brain space

% realign  (registration - align anat <-> epi)           
!@auto_tlrc -base MNI_avg152T1+tlrc -dxyz 3 -input anat+orig -ok_notice

% warp resting_ss using anat_ns in MNI space
!adwarp -apar anat+tlrc -dpar rest_ss+orig -dxyz 3 -prefix rest_warped
!3dAFNItoNIFTI rest_warped+tlrc rest_warped


                 
% ( 5 ) smoothing-------------------------------Gaussian blurring -------------
!3dmerge -1blur_fwhm 10.0 -doall -prefix rest_blr.nii rest_warped.nii



%%% --- resample rest_blr to match Smith3x3 templates
!3dresample -master ./RSN_HC_unthresholded_tmaps.nii -prefix master_cut.nii -inset rest_blr.nii





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - gift_ica on master_cut.nii.      Requires GIFT_ICA to be installed and on the matlab path
%                                    http://mialab.mrn.org/software/gift/index.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- make gift output directory
mkdir('gift_output');


% copy Input_data_subjects_1.m parameter file into subject directory. This
% file is provided by GIFT_ICA and can be edited manually to specify ICA
% parameters
%!cp Input_data_subjects_1.m ./

% run gift ica using icatb_batch_file_run
%icatb_batch_file_run('Input_data_subjects_1.m');      % uses file master_cut.nii

