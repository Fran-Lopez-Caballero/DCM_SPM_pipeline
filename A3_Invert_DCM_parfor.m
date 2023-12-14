%% Instructions

% 1) Define variables inside the parfor loop: EEG/MEG, NMDA/CMC, tag, etc.
% 2) Run one combination (EEG, CMC, pMMN, etc.) at a time

%% Define participant variable and load SPM
running_in = 'server'; % 'server' or 'local'

if strcmp(running_in,'server')
    addpath('/Private/User/spm12_v7771');
    addpath('/data/CNRL/data03/PEPP/User/DCM/Scripts_final_pipeline');
    subject_file_path = '/data/CNRL/data03/PEPP/User/DCM/Scripts_final_pipeline';
    outcome_path = '/data/CNRL/data03/PEPP/User/DCM/SPM_MEG_data';
    addpath('/data/CNRL/data03/PEPP/User/DCM/Scripts');
elseif strcmp(running_in,'local')
    addpath('C:/Users/user/private/path/Documents/MATLAB/spm12');
    addpath('C:/Private/User/DCM/Scripts_final_pipeline');
    subject_file_path = 'C:/Private/User/DCM/Scripts_final_pipeline';
    outcome_path = 'C:/Private/User/DCM/SPM_MEG_data';
    addpath('C:/Private/User/DCM/Scripts');
end

% participant_group = {'C','FE'};
subject_array_filename = 'subject_array_EEG_matched'; % subject_array_matched_baseline subject_array_EEG_matched
load([subject_file_path '/' subject_array_filename '.mat'],'subject_array');
participant = {subject_array{:,1}}; %#ok<*CCAT1,*USENS>
spm('defaults', 'EEG');

%% Parfor to invert DCMs

% Specify inside the loop which DCMs you are going to invert

% Only if R2 from that subject is currently under the criteria will DCM be re-inverted with new priors
parfor i = 1:length(participant)
    
    modality_data = 'EEG'; % Choose one at a time ('EEG' or 'MEG'
    MMN_types = 'dMMN'; % Choose one at a time (pMMN or dMMN)
    DCM_type = 'CMC'; % 'CMC' or 'NMDA'
    % CRITERIA_R2 = 85; % participants with fits under this R2 will be attempted to be re-inversed
    tag_DCM_solved = '_v3';
    
    % Load DCMs with new priors to invert
    output = parload([outcome_path '/DCM_structs/DCMs_' modality_data '_' MMN_types '_' DCM_type '.mat']);
    DCMs = output.DCMs;
    % Load associated participant index
    output = parload([outcome_path '/DCM_structs/Participant_index_' modality_data '_' MMN_types '_' DCM_type '.mat']);
    participant_index = output.participant_index;
    
    % Load DCM for this subject (MUST BE SAME LENGTH THAN PARTICIPANT!)
    DCM = DCMs{i,1};
           
   disp(['Inverting DCM for ' participant_index{i,1} ' ' modality_data ' ' DCM_type ' ' MMN_types]);

   % Invert DCM with new priors
   [DCM] = spm_dcm_fit(DCM); %#ok<*NASGU>
   DCM = DCM{1,1}; % For whatever reason it puts it into a cell
   parsave([outcome_path '/' participant_index{i,1} '_' modality_data '_DCM_' DCM_type '_' MMN_types '_solved' tag_DCM_solved '.mat'],DCM,'DCM');

end


disp(' ');      
disp('-------------------------');  
disp('DONE'); 
disp(datetime)
disp('-------------------------');     
disp(' '); 
