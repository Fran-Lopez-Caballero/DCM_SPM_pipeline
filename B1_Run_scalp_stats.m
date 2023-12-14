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
subject_array_filename = 'subject_array_EEG_matched'; % subject_array_EEG_matched subject_array_matched_baseline.mat
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
clear matlabbatch

%% Define Group order
groups = {'FE', 'C'};
conds = 'DDev'; % DDEV PDev
% conds = {'DDev'}; % DDEV PDev
count=0;
eval(['cd ' outcome_path]) % To ensure it saves them there
subs1=strings; pos_1 = 1;
subs2=strings; pos_2 = 1;
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
        if strcmp(modality_data{1},'EEG')
            if strcmp(subject_array{pos_subj,2},'FE')
                subs1(pos_1) = [participant{p} '_' modality_data{1} '_MMN']; %#ok<*AGROW>
                pos_1 = pos_1 + 1;
            elseif strcmp(subject_array{pos_subj,2},'C')
                subs2(pos_2) = [participant{p} '_' modality_data{1} '_MMN']; %#ok<*AGROW>
                pos_2 = pos_2 + 1;
            end
        elseif strcmp(modality_data{1},'MEG')
            if strcmp(subject_array{pos_subj,2},'FE')
                subs1(pos_1) = ['C' participant{p} '_' modality_data{1} '_MMN']; %#ok<*AGROW>
                pos_1 = pos_1 + 1;
            elseif strcmp(subject_array{pos_subj,2},'C')
                subs2(pos_2) = ['C' participant{p} '_' modality_data{1} '_MMN']; %#ok<*AGROW>
                pos_2 = pos_2 + 1;
            end
        end
    end
end

matlabbatch{1}.spm.stats.factorial_design.dir = {'GAVR_scalp_stats'};
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'Subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'Group';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'Tone Type';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;

count=0;
for i = 1:length(subs1)
    count=count+1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(count).scans = cellstr(append('SPM_images\',subs1{i},'\S_condition_',conds,'.nii'))';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(count).conds = [1 1; 1 2];
end
for i = 1:length(subs2)
    count=count+1;
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(count).scans = cellstr(append('SPM_images\',subs2{i},'\S_condition_',conds,'.nii'))';
    matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.fsubject(count).conds = [2 1; 2 2];
end

%matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums = [2; 3;];
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums = [2; 3;];
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{2}.fmain.fnum = 1;

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

spm_jobman('run',matlabbatch);
clear matlabbatch

matlabbatch{1}.spm.stats.fmri_est.spmmat = {'GAVR_scalp_stats/SPM.mat'};
matlabbatch{1}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical = 1;
spm_jobman('run',matlabbatch);
clear matlabbatch

%% Set Contrasts
nsubs = length(subs1)+length(subs2);
% Std vs Dev
matlabbatch{1, 1}.spm.stats.con.spmmat = {'GAVR_scalp_stats/SPM.mat'};
matlabbatch{1, 1}.spm.stats.con.consess{1,1}.fcon.name='Std vs Dev (all)';
matlabbatch{1, 1}.spm.stats.con.consess{1,1}.fcon.weights=[1 -1 1 -1 zeros(1,nsubs)/nsubs];
matlabbatch{1, 1}.spm.stats.con.consess{1,1}.fcon.sessrep = 'none';
% Std (FE-Sz) vs Std (Matched Control)
matlabbatch{1, 1}.spm.stats.con.consess{1,2}.fcon.name='Std (FE) vs Std (C)';
matlabbatch{1, 1}.spm.stats.con.consess{1,2}.fcon.weights=[1 0 -1 0 ones(1,length(subs1))/length(subs1) -ones(1, length(subs2))/length(subs2)];
matlabbatch{1, 1}.spm.stats.con.consess{1,2}.fcon.sessrep = 'none';
% Dev (FE-Sz) vs Dev (Matched Control)
matlabbatch{1, 1}.spm.stats.con.consess{1,3}.fcon.name='Dev (FE) vs Dev (C)';
matlabbatch{1, 1}.spm.stats.con.consess{1,3}.fcon.weights=[0 1 0 -1 ones(1,length(subs1))/length(subs1) -ones(1, length(subs2))/length(subs2)];
matlabbatch{1, 1}.spm.stats.con.consess{1,3}.fcon.sessrep = 'none';
spm_jobman('run',matlabbatch);
clear matlabbatch
