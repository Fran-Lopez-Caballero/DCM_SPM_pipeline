%% Define variables and start SPM

running_in = 'local'; % 'server' or 'local'

if strcmp(running_in,'server')
    rmpath('/private/path/brainstorm3_v20220706'); % Remove because of overlapping fieltrip functions
    addpath('/private/path/spm12_v7771');
    addpath('private/path/User/DCM/Scripts_final_pipeline');
    subject_file_path = 'private/path/User/DCM/Scripts_final_pipeline';
    load('private/path/User/DCM/Scripts_final_pipeline/Labels_EEG.mat');
    root_dir_analysis_path = 'private/path/analysis';
    brainstorm_anatomy_path = 'private/path/brainstorm_db/PEPP_FEManat_baseline/anat'; % For Baseline
    HCProc_folder = 'private/path/HCPproc';
    root_dir_bs = 'private/path/brainstorm_db/PEPP_MMN/data';
    outcome_path = 'private/path/User/DCM/SPM_MEG_data';
    addpath('private/path/User/DCM/Scripts');
elseif strcmp(running_in,'local')
    rmpath('C:/Users/private/path/Documents/MATLAB/brainstorm3'); % Remove because of overlapping fieltrip functions
    addpath('C:/Users/private/path/Documents/MATLAB/spm12');
    addpath('C:/private/User/DCM/Scripts_final_pipeline');
    subject_file_path = 'C:/private/User/DCM/Scripts_final_pipeline';
    load('C:/private/User/DCM/Scripts_final_pipeline/Labels_EEG.mat');
    root_dir_analysis_path = 'C:/private/analysis';
    brainstorm_anatomy_path = 'C:/private/brainstorm_db/PEPP_FEManat_baseline/anat'; % For Baseline
    HCProc_folder = 'C:/private/HCPproc';
    root_dir_bs = 'C:/private/brainstorm_db/PEPP_MMN/data';
    outcome_path = 'C:/private/User/DCM/SPM_MEG_data';
    addpath('C:/private/User/DCM/Scripts');
end

participant_group = {'C','FE'};
subject_array_filename = 'subject_array_EEG_matched'; % subject_array_matched_baseline subject_array_EEG_matched
load([subject_file_path '/' subject_array_filename '.mat']);
modality_data = {'EEG'}; % 'EEG', 'MEG', 'BIMODAL'
modify_triggers = 'YES'; % 'YES' or 'NO' Leave STI option or replace using vmrk file used in brainstorm
combine_planar = 'NO'; % 'YES' or 'NO' keep at NO so that further steps don't take too long
delete_previous_steps = 'YES'; % Deleting intermediate steps as they are created
delete_merged_files = 'YES'; % These have all trials without rejection, better to always keep
delete_average_files = 'YES'; % 'YES' or 'NO' MMN file has the same, so yes is ok, but since it takes a while to compute, just in case, no
crit_sweeps = 30; % minimum number of surviving sweeps to discard EEG, MEG or BIMODAL
crit_percent = 50; % minimum percentage of surviving sweeps to discard EEG, MEG or BIMODAL
fiducials_option = 'brainstorm'; % 'select' OR 'brainstorm'
EEG_sensor_combinations = {'EEG'}; % Only one
MEG_sensor_combinations = {'MEGandPLANAR','MEG','PLANAR'}; % In order based on Inverse_solution job numerical order for display
BIMODAL_sensor_combinations = {'EEGandMEGandPLANAR','EEGandMEG','EEGandPLANAR'}; % In order based on Inverse_solution job numerical order for display
compute_grand_average = 'YES'; % 'YES' or 'NO'
% Types of inverse solutions (explained in inverse solution section)
Inverse_solution_algorithm = {'IID','GS'}; % In order based on Inverse_solution job numerical order for display
type_of_average = 'standard'; % 'standard' or 'robust'

participant = {subject_array{:,1}};
spm('defaults', 'EEG');
condition_name = {'STD','PDEV','DDEV'};
condition_name_trigger = {'Standard','PDev','DDev'};

%% Convert to images

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,20},'ready_to_img')
        
        % Do amplitude thresholds separately for each modality
        for mod = 1:length(modality_data)
            
            % There will be no EEG file to average sweeps from
            if strcmp(modality_data{mod},'EEG') && (strcmp(subject_array{pos_subj,4},'no_EEG_chans') || strcmp(subject_array{pos_subj,5},'bad_EEG'))
                continue;
            end
            if strcmp(modality_data{mod},'MEG') && strcmp(subject_array{pos_subj,6},'bad_MEG')
                continue;
            end
            if strcmp(modality_data{mod},'BIMODAL') && (strcmp(subject_array{pos_subj,4},'no_EEG_chans') || strcmp(subject_array{pos_subj,5},'bad_EEG')...
                    || strcmp(subject_array{pos_subj,6},'bad_MEG') || strcmp(subject_array{pos_subj,7},'bad_BIMODAL'))
                continue;
            end
            
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},[modality_data{mod} '_MMN.mat']) & (...
                        startsWith({folders.name},[participant{p}])));
            if isempty(infolder)
                error(['No ' modality_data{mod} ' MMN files for ' participant{p}]);
            end

            % Should not be by definition, as it is 
            if length(infolder) > 1
                error(['More than one ' modality_data{mod} ' MMN file for ' participant{p}]);
            end

            eval(['cd ' outcome_path]) % To ensure it saves them there
            
            % Create nifti images
            for cond = 1:length(condition_name)                                
                matlabbatch{1}.spm.meeg.images.convert2images.D = cellstr([outcome_path '/' folders(infolder).name]);
                matlabbatch{1}.spm.meeg.images.convert2images.mode = 'scalp x time';
                matlabbatch{1}.spm.meeg.images.convert2images.conditions = {'Standard','PDev','DDev'};
                matlabbatch{1}.spm.meeg.images.convert2images.channels{1}.all = 'EEG';
                matlabbatch{1}.spm.meeg.images.convert2images.timewin = [-Inf Inf];
                matlabbatch{1}.spm.meeg.images.convert2images.freqwin = [-Inf Inf];
                matlabbatch{1}.spm.meeg.images.convert2images.prefix = '/SPM_images/';
                spm_jobman('run',matlabbatch);
                clear matlabbatch
            end
            
            
            for cond = 1:length(condition_name)
            folders = dir([outcome_path '/SPM_images/' participant{p} '_' modality_data{mod} '_MMN']);
            infolder = find(strcmp({folders.name},['condition_' condition_name_trigger{cond} '.nii']));
            if isempty(infolder)
                error(['No ' modality_data{mod} ' ' condition_name{cond} ' nii files for ' participant{p}]);
            end

            % Should not be by definition, as it is 
            if length(infolder) > 1
                error(['More than one ' modality_data{mod} ' ' condition_name{cond} ' nii file for ' participant{p}]);
            end
            
            % Smooth
            for cond = 1:length(condition_name)                                
                matlabbatch{1}.spm.spatial.smooth.data = cellstr([outcome_path '/SPM_images/' participant{p} '_' modality_data{mod} '_MMN/' folders(infolder).name]);
                matlabbatch{1}.spm.spatial.smooth.fwhm = [8 8 8];
                matlabbatch{1}.spm.spatial.smooth.dtype = 0;
                matlabbatch{1}.spm.spatial.smooth.im = 1;
                matlabbatch{1}.spm.spatial.smooth.prefix = 'S_';
                spm_jobman('run',matlabbatch);
                clear matlabbatch
            end
            end
            
            % No need to adjust names (resulting file is named the same)
            disp(' ');      
            disp('-------------------------');
            disp([participant{p} ' READY FOR STATS ' modality_data{mod}]);
            disp(datetime)
            disp(' ');                 
            
        end

    end   
end

%% Create matrix

eval(['cd ' outcome_path])
namesTemp = [];FE=[];C=[];
for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,20},'STATS')
        if strcmp(modality_data{1},'EEG') && (strcmp(subject_array{pos_subj,4},'no_EEG_chans') || strcmp(subject_array{pos_subj,5},'bad_EEG'))
            continue;
        end
        if strcmp(modality_data{1},'MEG') && strcmp(subject_array{pos_subj,6},'bad_MEG')
            continue;
        end
        if strcmp(modality_data{1},'BIMODAL') && (strcmp(subject_array{pos_subj,4},'no_EEG_chans') || strcmp(subject_array{pos_subj,5},'bad_EEG')...
                || strcmp(subject_array{pos_subj,6},'bad_MEG') || strcmp(subject_array{pos_subj,7},'bad_BIMODAL'))
            continue;
        end
        namesTemp{p,1} = [outcome_path '/' participant{p} '_' modality_data{1} '_MMN.mat']; %#ok<*SAGROW>
        if strcmp(subject_array{pos_subj,2},'FE')
            FE = [FE,p]; %#ok<*AGROW>
        elseif strcmp(subject_array{pos_subj,2},'C')
            C = [C,p];
        end
    end
end

%% Create Grand Means
clear matlabbatch

matlabbatch{1}.spm.meeg.averaging.grandmean.D = namesTemp(FE);
matlabbatch{1}.spm.meeg.averaging.grandmean.outfile = [modality_data{1} '_grandmean_FE'];
matlabbatch{1}.spm.meeg.averaging.grandmean.weighted = 0;
spm_jobman('run',matlabbatch);
clear matlabbatch

matlabbatch{1}.spm.meeg.averaging.grandmean.D = namesTemp(C);
matlabbatch{1}.spm.meeg.averaging.grandmean.outfile = [modality_data{1} '_grandmean_C'];
matlabbatch{1}.spm.meeg.averaging.grandmean.weighted = 0;
spm_jobman('run',matlabbatch);
clear matlabbatch
