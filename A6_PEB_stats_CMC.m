%% Statistics PEB DCM

%% Define variables and open SPM

clear;
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

% participant_group = {'C','FE'};
subject_array_filename = 'subject_array_EEG_matched'; % subject_array_matched_baseline subject_array_EEG_matched
load([subject_file_path '/' subject_array_filename '.mat']);
modality_data = 'EEG'; % Choose one at a time ('EEG' or 'MEG'
MMN_types = 'dMMN'; % Choose one at a time (pMMN or dMMN)
participant = {subject_array{:,1}};
tag_DCM_solved = '_v3';
spm('defaults', 'EEG');

%% Run PEB (Bayesian general linear model)

% Prepare input variables for PEB
GCM = {}; M = struct(); M.X = [];
for p = 1:length(participant)
    
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(modality_data,'EEG') && (strcmp(subject_array{pos_subj,4},'no_EEG_chans') || strcmp(subject_array{pos_subj,5},'bad_EEG'))
        continue;
    end
    if strcmp(modality_data,'MEG') && strcmp(subject_array{pos_subj,6},'bad_MEG')
        continue;
    end

    % Adjustement for 2379 not having pMMN MEG DCM
    if strcmp(modality_data,'MEG') && strcmp(MMN_types,'pMMN') && strcmp(subject_array{pos_subj,14},'no_pMMN_MEG_CMC_DCM')
        continue;
    end
    
    % Load DCM file for that subject
    try
        load([outcome_path '/' participant{p} '_' modality_data '_DCM_CMC_' MMN_types '_solved' tag_DCM_solved '.mat']);
    catch
        error(['No ' participant{p} '_' modality_data '_DCM_CMC_' MMN_types '_solved' tag_DCM_solved '.mat file found']);
    end
    
    M.X(p,1) = 1; % First column of all ones
    if contains(subject_array{pos_subj,2},'FE')
        M.X(p,2) = 1;
    elseif contains(subject_array{pos_subj,2},'C')
        M.X(p,2) = -1; 
    end
    GCM{p,1} = DCM;
end
field = {'all'};

% Eliminate empty cells (common issue with EEG/MEG alternative)
GCM = GCM(~cellfun('isempty',GCM));

% Do the same with M.X variable
M.X(~any(M.X,2),:) = [];  %rows

% Zach parameters (comment out to find original)
M.Q = 'all';
M.Xnames = {'Group','First-episode'};
M.Maxit=256; % Does not change results
field = {'A','B','M','N','G','T'}; % Changes results (less significant parameters when selecting only these)

% Run SPM function 
[PEB, DCM] = spm_dcm_peb(GCM,M,field);

% Peform "parameter prunning"
BMA = spm_dcm_peb_bmc(PEB);

% Inspect results with SPM GUI
spm_dcm_peb_review(BMA,GCM)
