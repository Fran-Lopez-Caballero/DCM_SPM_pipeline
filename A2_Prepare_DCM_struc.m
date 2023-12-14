%% Create DCMs with SPM analyzed data

% 1) Data must have been analyzed in SPM previously, with Analyze_in_SPM script
% 2) Inputs are files ending in MMN_DCM.mat and .dat
% 3) Subject array here must have ONLY the subjects you want DCM on

%% Define variables and start SPM

running_in = 'server'; % 'server' or 'local'

if strcmp(running_in,'server')
    % rmpath('/private/path/brainstorm3_v20220706'); % Remove because of overlapping fieltrip functions
    addpath('/private/path/spm12_v7771');
    addpath('~/matlab/brainstorm3_v20220706');
    addpath('private/path/User/DCM/Scripts/Analysis_in_SPM');
    subject_file_path = 'private/path/User/DCM/Scripts/Analysis_in_SPM';
    addpath('private/path/User/DCM/Scripts_final_pipeline');
    load('private/path/User/DCM/Scripts/Analysis_in_SPM/Labels_EEG.mat');
    root_dir_analysis_path = 'private/path/analysis';
    brainstorm_anatomy_path = 'private/path/brainstorm_db/PEPP_FEManat_baseline/anat'; % For Baseline
    HCProc_folder = 'private/path/HCPproc';
    root_dir_bs = 'private/path/brainstorm_db/PEPP_MMN/data';
    outcome_path = 'private/path/User/DCM/SPM_MEG_data';
    addpath('private/path/User/DCM/Scripts');
elseif strcmp(running_in,'local')
    rmpath('C:/Users/user/private/path/Documents/MATLAB/brainstorm3'); % Remove because of overlapping fieltrip functions
    addpath('C:/Users/user/private/path/Documents/MATLAB/spm12');
    addpath('private/path/User/DCM/Scripts/Analysis_in_SPM');
    subject_file_path = 'private/path/User/DCM/Scripts/Analysis_in_SPM';
    load('private/path/User/DCM/Scripts/Analysis_in_SPM/Labels_EEG.mat');
    root_dir_analysis_path = 'private/path/analysis';
    brainstorm_anatomy_path = 'private/path/brainstorm_db/PEPP_FEManat_baseline/anat'; % For Baseline
    HCProc_folder = 'private/path/HCPproc';
    root_dir_bs = 'private/path/brainstorm_db/PEPP_MMN/data';
    outcome_path = 'private/path/User/DCM/SPM_MEG_data';
    addpath('private/path/User/DCM/Scripts');
end

participant_group = {'C','FE'};
subject_array_filename = 'subject_array_EEG_matched';
load([subject_file_path '/' subject_array_filename '.mat']);
modality_data = {'EEG'}; % 'EEG', 'MEG', 'BIMODAL'
modify_triggers = 'YES'; % 'YES' or 'NO' Leave STI option or replace using vmrk file used in brainstorm
combine_planar = 'NO'; % 'YES' or 'NO' keep at NO so that further steps don't take too long
delete_previous_steps = 'NO'; % Deleting intermediate steps as they are created
delete_merged_files = 'YES'; % These have all trials without rejection, better to always keep
delete_average_files = 'YES'; % 'YES' or 'NO' MMN file has the same, so yes is ok, but since it takes a while to compute, just in case, no
crit_sweeps = 30; % minimum number of surviving sweeps to discard EEG, MEG or BIMODAL
crit_percent = 50; % minimum percentage of surviving sweeps to discard EEG, MEG or BIMODAL
fiducials_option = 'brainstorm'; % 'select' OR 'brainstorm'
EEG_sensor_combinations = {'EEG'}; % Only one
MEG_sensor_combinations = {'MEGandPLANAR','MEG','PLANAR'}; % In order based on Inverse_solution job numerical order for display
BIMODAL_sensor_combinations = {'EEGandMEGandPLANAR','EEGandMEG','EEGandPLANAR'}; % In order based on Inverse_solution job numerical order for display
% Types of inverse solutions (explained in inverse solution section)
Inverse_solution_algorithm = {'IID','GS'}; % In order based on Inverse_solution job numerical order for display
type_of_average = 'standard'; % 'standard' or 'robust'

MMN_types = {'pMMN','dMMN'};
MMN_combinations_DCM = {[1,2],[1,3]}; % 1 = STD, 2 = PDev, 3 = DDev
rois = {'left_A1','right_A1','left_STG','right_STG','left_IFG','right_IFG'};
mni_left_A1 = [-42, -22,7];
mni_right_A1 =  [46, -14, 8];
mni_left_STG = [-61, -32, 8];
mni_right_STG = [59, -25, 8];
mni_left_IFG = [-46, 20, 8];
mni_right_IFG = [46, 20, 8];
which_MEG_sensors = 'MEGMAG'; % DCM can only be run with one type ('MEGMAG' or 'MEGPLANAR')
participant = {subject_array{:,1}};
tag_DCM_solved = 'solved_v3';
spm('defaults', 'EEG');

%% Create DCM structures (CMC)

if strcmp(running_in,'local')
    error('This section cannot be run in local (personalized functions in server');
end

% VERY IMPORTANT: only one modality can be modeled at a time (MEG PLANAR, MEG or EEG)
[model_space] = gen_model_space();

for mod = 1:length(modality_data)
for mmnty = 1:length(MMN_types)
    DCMs = {};
    participant_index = {};
    for p = 1:length(participant)
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},[modality_data{mod} '_MMN_DCM.mat']) & (...
                    startsWith({folders.name},[participant{p}])));     
        if isempty(infolder)
            error(['No ' modality_data{mod} ' MMN DCM files for ' participant{p}]);
        end

        % Should not be by definition, as it is 
        if length(infolder) > 1
            error(['More than one ' modality_data{mod} ' MMN DCM file for ' participant{p}]);
        end

        % MMN DCM file
        DCM.xY.Dfile = ([outcome_path '/' folders(infolder).name]);
        MEEG = load(DCM.xY.Dfile);

        % Remove either PLANAR or MAG sensors so that it does not ask
        % (USELESS, WE NEED TO CROP BEFORE)
        % if strcmp(modality_data{mod},'MEG')
        %    % Change sensor labels so that it does not take that into account
        %    for ch = 1:length(MEEG.D.channels)
        %      if ~strcmp(MEEG.D.channels(ch).type,which_MEG_sensors)
        %           MEEG.D.channels(ch).type = 'DO_NOT_USE';
        %       end
        %   end                         
        % end
        D = MEEG.D;

        % COMMENTED SECTION AS WE ARE USING INDIVIDUAL MRIs already in the file    
        %    % Fix location of canonical brain in MEEG file
        %    %--------------------------------------------------------------------------
        %     D.path = D_file_path;
        %     D.data.fname = [D_file_path '\' participant{p} '_SPM.dat'];
        %     D.other.inv{1,1}.forward.vol = [canonical_path 'single_subj_T1_EEG_BEM.mat'];
        %     D.other.inv{1,1}.mesh.sMRI = [canonical_path 'single_subj_T1.nii'];
        %     D.other.inv{1,1}.mesh.tess_ctx = [canonical_path 'cortex_8196.surf.gii'];
        %     D.other.inv{1,1}.mesh.tess_scalp = [canonical_path 'scalp_2562.surf.gii'];
        %     D.other.inv{1,1}.mesh.tess_oskull = [canonical_path 'oskull_2562.surf.gii'];
        %     D.other.inv{1,1}.mesh.tess_iskull = [canonical_path 'iskull_2562.surf.gii'];
        %     save(DCM.xY.Dfile, 'D');

        % Parameters and options used for setting up model
        DCM.options.analysis = 'ERP'; % analyze evoked responses
        DCM.options.model    = 'CMC'; % CMC model
        DCM.options.spatial  = 'ECD'; % spatial model
        DCM.options.Tdcm(1)  = -50;     % start of peri-stimulus time to be modelled
        DCM.options.Tdcm(2)  = 300;   % end of peri-stimulus time to be modelled
        DCM.options.Nmodes   = 8;     % nr of modes for data selection
        DCM.options.h        = 1;     % nr of DCT components
        DCM.options.onset    = 60;    % selection of onset (prior mean) deafult = 75 USED 110
        DCM.options.D        = 1;     % downsampling
        DCM.options.location = 0;     % optimising source location default = 1
        DCM.options.han      = 1;     % applying hanning window
        DCM.options.trials   = MMN_combinations_DCM{mmnty};     % index of all ERPs within file  %

        % Data and spatial model
        DCM  = spm_dcm_erp_data(DCM);

        % Location priors for dipoles
        % Locations of MMN sources in template are based on Juanita Todd Schz Bulleting paper (2022):
        % "The boundary element model in Montreal Neurological Institute space was used as per the default SPM12 template. 
        % The modeled network assumed six cortical sources identified as integral to auditory inference40 
        % (figure 1A with MNI locations as follows: 
        % left A1 [-42, -22,7], right A1 [46, -14, 8], 
        % left STG [?61, ?32, 8], right STG [59, ?25, 8], 
        % left IFG [?46, 20, 8], right IFG [46, 20, 8])."

        anat_folder_bs = dir(brainstorm_anatomy_path);
        infolder_anat_bs = find(strcmp({anat_folder_bs.name},participant{p}));
        if isempty(infolder_anat_bs)
            error(['No bs anat folder for ' participant{p}]);
        end
        if length(infolder_anat_bs) > 1
            error(['More than one anat folder for ' participant{p}]);
        end
        load([brainstorm_anatomy_path '/' anat_folder_bs(infolder_anat_bs).name '/brainstormsubject.mat'],'Anatomy');
        % Retrieve the mni coordinates
        sMri = load([brainstorm_anatomy_path '/' anat_folder_bs(infolder_anat_bs).name '/' Anatomy(7:end)]);

        % Create n by 3 matrix for dipoles in native MRI space
        Dipole_matrix = [];
        for roi = 1:length(rois)
            eval(['current_coordinates = mni_' rois{roi} ';'])
            Dipole_matrix(:,roi) = cs_convert(sMri, 'mni', 'world', current_coordinates ./ 1000) .* 1000;
        end

        DCM.Lpos  = Dipole_matrix;
        DCM.Sname = {'left A1', 'right A1', 'left STG', 'right STG', 'left IFG', 'right IFG'};
        Nareas    = size(DCM.Lpos,2);

        % Spatial model
        DCM = spm_dcm_erp_dipfit(DCM);

        % Specify connectivity model
        DCM.A{1} = zeros(Nareas, Nareas);   % forward connections
        DCM.A{1}(3,1) = 1;
        DCM.A{1}(4,2) = 1;
        DCM.A{1}(5,3) = 1;
        DCM.A{1}(6,4) = 1;

        DCM.A{2} = zeros(Nareas,Nareas);    % backward connections
        DCM.A{2}(1,3) = 1;
        DCM.A{2}(2,4) = 1;
        DCM.A{2}(3,5) = 1;
        DCM.A{2}(4,6) = 1;

        DCM.A{3} = zeros(Nareas,Nareas);    % modulatory connections
        DCM.A{3}(1,1) = 1;
        DCM.A{3}(2,2) = 1;
        DCM.A{3}(3,3) = 1;
        DCM.A{3}(4,4) = 1;
        DCM.A{3}(5,5) = 1;
        DCM.A{3}(6,6) = 1;
        DCM.B{1} = zeros(Nareas,Nareas);
        DCM.B{1} = model_space{end}.matrix;    % model specification --  fully connected model inverted
        DCM.C = [1; 1; 0; 0; 0; 0];            % input
        % Between trial effects
        DCM.xU.X = [0 1]';
        DCM.xU.name{1} = {'Deviance'};

        % Do priors as well?
        [pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
        DCM.M.pC = pC;
        DCM.M.pE = pE;

        if strcmp(modality_data{mod},'EEG')
            DCM.M.Nmax = 64;
        elseif strcmp(modality_data{mod},'MEG')
            DCM.M.Nmax = 102; % ??? (MEG PLANAR/MEG)
        end

        % Define name (apparently it need to be different as it is called in the future)
        DCM.name = strcat([outcome_path '/' participant{p} '_' modality_data{mod} '_DCM_CMC_' MMN_types{mmnty} '_' tag_DCM_solved '.mat']);

        disp(' ');      
        disp('-------------------------');
        disp(['Creating DCM CMC for ' participant{p} '(' modality_data{mod} ' ' MMN_types{mmnty} ')...']);
        disp(datetime)
        disp(' '); 

        DCMs{p,1} = DCM; %#ok<*NODEF,*AGROW> 
        participant_index{p,1} = participant{p}; %#ok<*SAGROW> 
        clear('DCM');
    end
    DCMs = DCMs(~cellfun('isempty',DCMs)); % In case there is any empty
    participant_index = participant_index(~cellfun('isempty',participant_index));
    if length(participant_index) ~= length(DCMs)
        error('participant list does not match DCM list');
    end
    % Save variables
    save([outcome_path '/DCM_structs/DCMs_' modality_data{mod} '_' MMN_types{mmnty} '_CMC.mat'],'DCMs'); 
    save([outcome_path '/DCM_structs/Participant_index_' modality_data{mod} '_' MMN_types{mmnty} '_CMC.mat'],'participant_index');
end
end

disp('DONE creating DCMs for CMC');
 
%% Create DCM structures(NMDA)

if strcmp(running_in,'local')
    error('This section cannot be run in local (personalized functions in server');
end

% VERY IMPORTANT: only one modality can be modeled at a time (MEG PLANAR, MEG or EEG)
[model_space] = gen_model_space();

for mod = 1:length(modality_data)
for mmnty = 1:length(MMN_types)
    DCMs = {};
    participant_index = {};
    for p = 1:length(participant)
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},[modality_data{mod} '_MMN_DCM.mat']) & (...
                    startsWith({folders.name},[participant{p}])));
        if isempty(infolder)
            error(['No ' modality_data{mod} ' MMN DCM files for ' participant{p}]);
        end

        % Should not be by definition, as it is 
        if length(infolder) > 1
            error(['More than one ' modality_data{mod} ' MMN DCM file for ' participant{p}]);
        end

        % MMN DCM file
        DCM.xY.Dfile = ([outcome_path '/' folders(infolder).name]);
        MEEG                                = load(DCM.xY.Dfile);
        D                                   = MEEG.D;

        %    COMMENTED SECTION AS WE ARE USING INDIVIDUAL MRIs already in the file    
        %    % Fix location of canonical brain in MEEG file
        %    %--------------------------------------------------------------------------
        %     D.path = D_file_path;
        %     D.data.fname = [D_file_path '\' participant{p} '_SPM.dat'];
        %     D.other.inv{1,1}.forward.vol = [canonical_path 'single_subj_T1_EEG_BEM.mat'];
        %     D.other.inv{1,1}.mesh.sMRI = [canonical_path 'single_subj_T1.nii'];
        %     D.other.inv{1,1}.mesh.tess_ctx = [canonical_path 'cortex_8196.surf.gii'];
        %     D.other.inv{1,1}.mesh.tess_scalp = [canonical_path 'scalp_2562.surf.gii'];
        %     D.other.inv{1,1}.mesh.tess_oskull = [canonical_path 'oskull_2562.surf.gii'];
        %     D.other.inv{1,1}.mesh.tess_iskull = [canonical_path 'iskull_2562.surf.gii'];
        %     save(DCM.xY.Dfile, 'D');

        % Parameters and options used for setting up model
        DCM.options.analysis = 'ERP'; % analyze evoked responses
        DCM.options.model    = 'CMM_NMDA'; % CMC model
        DCM.options.spatial  = 'ECD'; % spatial model
        DCM.options.Tdcm(1)  = -50;     % start of peri-stimulus time to be modelled
        DCM.options.Tdcm(2)  = 300;   % end of peri-stimulus time to be modelled
        DCM.options.Nmodes   = 8;     % nr of modes for data selection
        DCM.options.h        = 1;     % nr of DCT components
        DCM.options.onset    = 60;    % selection of onset (prior mean) deafult = 75 USED 110
        DCM.options.D        = 1;     % downsampling
        DCM.options.location = 0;     % optimising source location default = 1
        DCM.options.han      = 1;     % applying hanning window
        DCM.options.trials   = MMN_combinations_DCM{mmnty};     % index of all ERPs within file  %

        % Data and spatial model
        DCM  = spm_dcm_erp_data(DCM);

        % Location priors for dipoles    
        % Locations of MMN sources in template are based on Juanita Todd Schz Bulleting paper (2022):
        % "The boundary element model in Montreal Neurological Institute space was used as per the default SPM12 template. 
        % The modeled network assumed six cortical sources identified as integral to auditory inference40 
        % (figure 1A with MNI locations as follows: 
        % left A1 [-42, -22,7], right A1 [46, -14, 8], 
        % left STG [?61, ?32, 8], right STG [59, ?25, 8], 
        % left IFG [?46, 20, 8], right IFG [46, 20, 8])."

        anat_folder_bs = dir(brainstorm_anatomy_path);
        infolder_anat_bs = find(strcmp({anat_folder_bs.name},participant{p}));
        if isempty(infolder_anat_bs)
            error(['No bs anat folder for ' participant{p}]);
        end
        if length(infolder_anat_bs) > 1
            error(['More than one anat folder for ' participant{p}]);
        end
        load([brainstorm_anatomy_path '/' anat_folder_bs(infolder_anat_bs).name '/brainstormsubject.mat'],'Anatomy');
        % Retrieve the mni coordinates
        sMri = load([brainstorm_anatomy_path '/' anat_folder_bs(infolder_anat_bs).name '/' Anatomy(7:end)]);

        % Create n by 3 matrix for dipoles in native MRI space
        Dipole_matrix = [];
        for roi = 1:length(rois)
            eval(['current_coordinates = mni_' rois{roi} ';'])
            Dipole_matrix(:,roi) = cs_convert(sMri, 'mni', 'world', current_coordinates ./ 1000) .* 1000;
        end

        DCM.Lpos  = Dipole_matrix;
        DCM.Sname = {'left A1', 'right A1', 'left STG', 'right STG', 'left IFG', 'right IFG'};
        Nareas    = size(DCM.Lpos,2);

        % Spatial model
        DCM = spm_dcm_erp_dipfit(DCM);

        % Specify connectivity model
        DCM.A{1} = zeros(Nareas, Nareas);   % forward connections
        DCM.A{1}(3,1) = 1;
        DCM.A{1}(4,2) = 1;
        DCM.A{1}(5,3) = 1;
        DCM.A{1}(6,4) = 1;

        DCM.A{2} = zeros(Nareas,Nareas);    % backward connections
        DCM.A{2}(1,3) = 1;
        DCM.A{2}(2,4) = 1;
        DCM.A{2}(3,5) = 1;
        DCM.A{2}(4,6) = 1;

        DCM.A{3} = zeros(Nareas,Nareas);    % modulatory connections
        % DCM.A{3}(1,1) = 1; % Commented since this is NMDA
        % DCM.A{3}(2,2) = 1;
        % DCM.A{3}(3,3) = 1;
        % DCM.A{3}(4,4) = 1;
        % DCM.A{3}(5,5) = 1;
        % DCM.A{3}(6,6) = 1;
        DCM.B{1} = zeros(Nareas,Nareas);
        DCM.B{1}(3,1) = 1;
        DCM.B{1}(4,2) = 1;
        DCM.B{1}(5,3) = 1;
        DCM.B{1}(6,4) = 1;
        DCM.B{1}(1,3) = 1;
        DCM.B{1}(2,4) = 1;
        DCM.B{1}(3,5) = 1;
        DCM.B{1}(4,6) = 1;
        % DCM.B{1} = model_space{end}.matrix;           % model specification --  fully connected model inverted
        DCM.C = [1; 1; 0; 0; 0; 0]; % input
        % Between trial effects
        DCM.xU.X = [0 1]';
        DCM.xU.name{1} = {'Deviance'};

        % Do priors as well?
        [pE,pC]  = spm_dcm_neural_priors(DCM.A,DCM.B,DCM.C,DCM.options.model);
        DCM.M.pC = pC;
        DCM.M.pE = pE;

        if strcmp(modality_data{mod},'EEG')
            DCM.M.Nmax = 64;
        elseif strcmp(modality_data{mod},'MEG')
            DCM.M.Nmax = 102; % ??? (MEG PLANAR/MEG)
        end

        % Define name
        DCM.name = strcat([outcome_path '/' participant{p} '_' modality_data{mod} '_DCM_NMDA_' MMN_types{mmnty} '_' tag_DCM_solved '.mat']);

        disp(' ');      
        disp('-------------------------');
        disp(['Creating DCM NMDA for ' participant{p} ' (' modality_data{mod} ' ' MMN_types{mmnty} ')...']);
        disp(datetime)
        disp(' '); 

        DCMs{p,1} = DCM; %#ok<*NODEF,*AGROW> 
        participant_index{p,1} = participant{p}; %#ok<*SAGROW>         
        clear('DCM');
    end
    DCMs = DCMs(~cellfun('isempty',DCMs)); % In case there is any empty
    participant_index = participant_index(~cellfun('isempty',participant_index));

    if length(participant_index) ~= length(DCMs)
        error('participant list does not match DCM list');
    end
    % Save variables
    save([outcome_path '/DCM_structs/DCMs_' modality_data{mod} '_' MMN_types{mmnty} '_NMDA.mat'],'DCMs'); 
    save([outcome_path '/DCM_structs/Participant_index_' modality_data{mod} '_' MMN_types{mmnty} '_NMDA.mat'],'participant_index');
end  
end


 disp('DONE creating DCMs for NMDA');
