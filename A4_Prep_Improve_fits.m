%% Define variables
running_in = 'local'; % 'server' or 'local'

if strcmp(running_in,'server')
    addpath('/Private/User/spm12_v7771');
    addpath('/private/path/private/path/Fran/DCM/Scripts_final_pipeline');
    subject_file_path = '/private/path/private/path/Fran/DCM/Scripts_final_pipeline';
    outcome_path = '/private/path/private/path/Fran/DCM/SPM_MEG_data';
    addpath('/private/path/private/path/Fran/DCM/Scripts');
elseif strcmp(running_in,'local')
    addpath('C:/Users/USER/private/Documents/MATLAB/spm12');
    addpath('C:/private/path/Fran/DCM/Scripts_final_pipeline');
    subject_file_path = 'C:/private/path/Fran/DCM/Scripts_final_pipeline';
    outcome_path = 'C:/private/path/Fran/DCM/SPM_MEG_data';
    addpath('C:/private/path/Fran/DCM/Scripts');
end

% participant_group = {'C','FE'};
subject_array_filename = 'subject_array_EEG_matched'; % subject_array_matched_baseline subject_array_EEG_matched
load([subject_file_path '/' subject_array_filename '.mat'],'subject_array');
modality_data = 'EEG'; % Choose one at a time ('EEG' or 'MEG'
MMN_types = 'dMMN'; % Choose one at a time (pMMN or dMMN)
DCM_type = 'CMC'; % 'CMC' or 'NMDA'
% R2=100*Predicted response/(Predicted response + Residual response). Response being PCA first eight components
participant = {subject_array{:,1}}; %#ok<*CCAT1,*USENS>
tag_DCM_solved = '_v3';
spm('defaults', 'EEG');

%% Define DCM and participant array

DCMs = {};
participant_index = {};
% Load solved 
for p = 1:length(participant)
    
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(modality_data,'EEG') && (strcmp(subject_array{pos_subj,4},'no_EEG_chans') || strcmp(subject_array{pos_subj,5},'bad_EEG'))
        continue;
    end
    if strcmp(modality_data,'MEG') && strcmp(subject_array{pos_subj,6},'bad_MEG')
        continue;
    end
    
    if strcmp(modality_data,'MEG') && strcmp(MMN_types,'pMMN') && strcmp(subject_array{pos_subj,14},['no_pMMN_MEG_' DCM_type '_DCM'])
        continue;
    end

    % Load DCM file for that subject
    try
       load([outcome_path '/' participant{p} '_' modality_data '_DCM_' DCM_type '_' MMN_types '_solved' tag_DCM_solved '.mat'],'DCM');
    catch
        error(['No ' participant{p} '_' modality_data '_DCM_' DCM_type '_' MMN_types '_solved' tag_DCM_solved '.mat file found']);
    end

    DCMs{p,1} = DCM; %#ok<*NODEF,*AGROW>
    participant_index{p,1} = participant{p}; %#ok<*SAGROW>
end

% Eliminate empty cells (common issue with EEG/MEG alternative)
DCMs = DCMs(~cellfun('isempty',DCMs));
participant_index = participant_index(~cellfun('isempty',participant_index));

if length(participant_index) ~= length(DCMs)
    error('participant list does not match DCM list');
end

%% Compute R2 of each DCM

R2 = zeros(1,length(DCMs));
for i = 1:length(DCMs)
    RSS=sum(sum([DCMs{i}.R{1}.^2 DCMs{i}.R{2}.^2]));
    PSS=sum(sum([(DCMs{i}.H{1}+DCMs{i}.R{1}).^2 (DCMs{i}.H{2}+DCMs{i}.R{2}).^2]));
    R2(i)=100*PSS/(PSS + RSS);
    participant_index{i,2} = R2(i);
end

%% Update priors based on PEB

M   = struct();
M.Q = 'all';
M.X(:,1) = ones(length(DCMs),1);
field = {'all'};
% field = {'A{1}'};

[~,DCMs] = spm_dcm_peb(DCMs,M,field);
% [~,DCMs2] = spm_dcm_peb(DCMs,M,field);
%% Save structures for next script in pipeline

% Save all DCMs with new priors (regardless of whether re-inverted or not) for consistency
save([outcome_path '/Updated_priors/DCMs_' modality_data '_DCM_' DCM_type '_' MMN_types '.mat'],'DCMs');
% Also, save participant-index variable to load it in parloop
save([outcome_path '/Updated_priors/Participant_index_' modality_data '_DCM_' DCM_type '_' MMN_types '.mat'],'participant_index');
