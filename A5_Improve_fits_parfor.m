%% Instructions

% 1) Define variables inside the parfor loop: EEG/MEG, NMDA/CMC, tag, etc.
% 2) Run one combination (EEG, CMC, pMMN, etc.) at a time

%% Define participant variable and load SPM
running_in = 'server'; % 'server' or 'local'

if strcmp(running_in,'server')
    addpath('/Private/User/spm12_v7771');
    addpath('/private/path/private/path/Fran/DCM/Scripts_final_pipeline');
    subject_file_path = '/private/path/private/path/Fran/DCM/Scripts_final_pipeline';
    outcome_path = '/private/path/private/path/Fran/DCM/SPM_MEG_data';
    addpath('/private/path/private/path/Fran/DCM/Scripts');
    addpath('/private/path/private/DCM/Scripts_final_pipeline');
elseif strcmp(running_in,'local')
    addpath('C:/Users/USER/private/Documents/MATLAB/spm12');
    addpath('C:/private/path/Fran/DCM/Scripts_final_pipeline');
    subject_file_path = 'C:/private/path/Fran/DCM/Scripts_final_pipeline';
    outcome_path = 'C:/private/path/Fran/DCM/SPM_MEG_data';
    addpath('C:/private/path/Fran/DCM/Scripts');
    addpath('/private/path/private/Scripts_final_pipeline');
end

% participant_group = {'C','FE'};
subject_array_filename = 'subject_array_EEG_matched'; % subject_array_matched_baseline subject_array_EEG_matched
load([subject_file_path '/' subject_array_filename '.mat'],'subject_array');
participant = {subject_array{:,1}}; %#ok<*CCAT1,*USENS>
cd '/data/CNRL/data03/PEPP/Fran/DCM/SPM_MEG_data';
spm('defaults', 'EEG');

%% Parfor to invert DCMs

% Only if R2 from that subject is currently under the criteria will DCM be re-inverted with new priors
parfor i = 1:length(participant)
    
    modality_data = 'EEG'; % Choose one at a time ('EEG' or 'MEG'
    MMN_types = 'dMMN'; % Choose one at a time (pMMN or dMMN)
    DCM_type = 'CMC'; % 'CMC' or 'NMDA'
    CRITERIA_R2 = 85; % participants with fits under this R2 will be attempted to be re-inversed
    tag_DCM_solved = '_v3';
    outcome_path = '/private/path/DCM/SPM_MEG_data';
    
    % Load DCMs with new priors to invert
    output = parload([outcome_path '/Updated_priors/DCMs_' modality_data '_DCM_' DCM_type '_' MMN_types '.mat']);
    DCMs = output.DCMs;
    % Load associated participant index
    output = parload([outcome_path '/Updated_priors/Participant_index_' modality_data '_DCM_' DCM_type '_' MMN_types '.mat']);
    participant_index = output.participant_index;
    
    % Load DCM for this subject (MUST BE SAME LENGTH THAN PARTICIPANT!)
    DCM = DCMs{i,1};
    
    % Save with new priors for consistency in naming (regardless of whether re-inverted or not)
    parsave([outcome_path '/' participant_index{i,1} '_' modality_data '_DCM_' DCM_type '_' MMN_types '_solved' tag_DCM_solved '.mat'],DCM,'DCM');
    
    if participant_index{i,2} < CRITERIA_R2  %#ok<*PFBNS>
       
       disp(['Inverting DCM for ' participant_index{i,1} ' ' modality_data ' ' DCM_type ' ' MMN_types]);
       
       % Invert DCM with new priors
       [DCM] = spm_dcm_fit(DCM); %#ok<*NASGU>
       DCM = DCM{1,1}; % For whatever reason it puts it into a cell
       
       % Calculate R2 again: if it's worst than before, don't save
       RSS=sum(sum([DCM.R{1}.^2 DCM.R{2}.^2]));
       PSS=sum(sum([(DCM.H{1}+DCM.R{1}).^2 (DCM.H{2}+DCM.R{2}).^2]));
       R2_new=100*PSS/(PSS + RSS);
       if R2_new > participant_index{i,2}
           % Overwrite solved 2 with better fits
           parsave([outcome_path '/' participant_index{i,1} '_' modality_data '_DCM_' DCM_type '_' MMN_types '_solved' tag_DCM_solved '.mat'],DCM,'DCM');
       end
    end
end


disp(' ');      
disp('-------------------------');  
disp('DONE'); 
disp(datetime)
disp('-------------------------');     
disp(' '); 
