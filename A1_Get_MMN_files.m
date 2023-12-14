%% Prepare MMN files to run DCM with SPM analyzed data

% * You must change subject array third column to 'needs_cropch' first

% * To run this script, S0_Analyze_in_SPM.m MUST have been ran first,
% which analyzed all data in SPM including sources using MNE. The files
% that we need from that script are the ones ending in cleaned.mat/dat,
% which have already all the preprocessing up to the sources

% * The first three sections (crop, average across conditions and obtain MMN) 
% help us separate files from the ones with same names obtaiend in the
% Compute_sources_in_SPM.m script. From there, we can run dipole fitting

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
spm('defaults', 'EEG');

%% Previous step 1) Crop channels that won't be used in each modality

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_cropch')
        
        % Do amplitude thresholds separately for each modality
        for mod = 1:length(modality_data)    
            
            % There will be no file to sort conditions with, or no necessity to do this as the file won't be used
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
            infolder = find(endsWith({folders.name},[modality_data{mod} '_cleaned.mat']) & (...
                        startsWith({folders.name},[participant{p}])));
            if isempty(infolder)
                error(['No ' modality_data{mod} ' cleaned files for ' participant{p}]);
            end

            % Should not be by definition, as it is 
            if length(infolder) > 1
                error(['More than one ' modality_data{mod} ' cleaned file for ' participant{p}]);
            end
            
            jobfile = {[subject_file_path '/Crop_' modality_data{mod} '_channels_job.m']};
            jobs = repmat(jobfile, 1, 1);
            inputs = cell(0, 1);
            inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]);
            eval(['cd ' outcome_path]) % To ensure it saves them there
            spm('defaults', 'EEG');
            spm_jobman('run', jobs, inputs{:});


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust names
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_cleaned.mat') & (...
                        startsWith({folders.name},['C' participant{p} '_' modality_data{mod}])));
            % Rename each run file
            for in = 1:length(infolder)
                line = infolder(in);
                load([outcome_path '/' folders(line).name]);
                D.data.fname = [outcome_path '/' participant{p} '_' modality_data{mod} '_cropped.dat'];
                D.fname = [participant{p} '_' modality_data{mod} '_cropped.mat'];

                % Adjust name of .dat file with data as well
                datfilename = [outcome_path '/C' participant{p} '_' modality_data{mod} '_cleaned.dat'];  %#ok<*SAGROW> 
                % If it was already renamed there will be not imported_dat
                if exist(datfilename,'file')
                    datnewname = [outcome_path '/' participant{p} '_' modality_data{mod} '_cropped.dat'];  %#ok<*SAGROW>  
                    movefile(datfilename,datnewname);
                end

                % Now, save .mat
                save([outcome_path '/' participant{p} '_' modality_data{mod} '_cropped'],'D');
                % Delete original file
                delete([outcome_path '/C' participant{p} '_' modality_data{mod} '_cleaned.mat']);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Here never delete the cleaned file
            
            disp(' ');      
            disp('-------------------------');
            disp([participant{p} ' CROPED ' modality_data{mod}]);
            disp(datetime)
            disp(' '); 
        end

        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_average';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 

    end   
end

%% Previous step 2) average across conditions (DEV, STD, MMN)

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_average')
        
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
            infolder = find(endsWith({folders.name},[modality_data{mod} '_cropped.mat']) & (...
                        startsWith({folders.name},[participant{p}])));
            if isempty(infolder)
                error(['No ' modality_data{mod} ' cropped files for ' participant{p}]);
            end

            % Should not be by definition, as it is 
            if length(infolder) > 1
                error(['More than one ' modality_data{mod} ' cropped file for ' participant{p}]);
            end

            jobfile = {[subject_file_path '/Average_' type_of_average '_job.m']};
            jobs = repmat(jobfile, 1, 1);
            inputs = cell(0, 1);
            inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]);
            eval(['cd ' outcome_path]) % To ensure it saves them there
            spm('defaults', 'EEG');
            spm_jobman('run', jobs, inputs{:});


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust names
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},[modality_data{mod} '_cropped.mat']) & (...
                        startsWith({folders.name},['M' participant{p}])));
            if length(infolder) > 1
                error(['More than one averaged file for ' participant{p} ' ' modality_data{mod}]);
            end      

            % Should only be one, but just in case
            load([outcome_path '/' folders(infolder).name]);
            D.data.fname = [outcome_path '/' participant{p} '_' modality_data{mod} '_averaged.dat'];
            D.fname = [participant{p} '_' modality_data{mod} '_averaged.mat'];

            % Adjust name of .dat file with data as well
            datfilename = [outcome_path '/M' participant{p} '_' modality_data{mod} '_cropped.dat'];  %#ok<*SAGROW> 
            % If it was already renamed there will be not imported_dat
            if exist(datfilename,'file')
                datnewname = [outcome_path '/' participant{p} '_' modality_data{mod} '_averaged.dat'];  %#ok<*SAGROW>  
                movefile(datfilename,datnewname);
            end

            % Also, delete M avergae file
            delete([outcome_path '/M' participant{p} '_' modality_data{mod} '_cropped.mat']);
            
            % Now, save .mat (will overwrite previous .mat)
            save([outcome_path '/' participant{p} '_' modality_data{mod} '_averaged'],'D');                
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Averaged file will already be cropped, so:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Delete previous step
            if strcmp(delete_previous_steps,'YES')
                delete([outcome_path '/' participant{p} '_' modality_data{mod} '_cropped.dat']);
                delete([outcome_path '/' participant{p} '_' modality_data{mod} '_cropped.mat']);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % No need to adjust names (resulting file is named the same)
            disp(' ');      
            disp('-------------------------');
            disp([participant{p} ' AVERAGED ' modality_data{mod}]);
            disp(datetime)
            disp(' ');                 
            
        end

        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_MMN';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
        
        % Do not delete previous steps (cleaned)

    end   
end

%% Obtain Mismatch but call it differently than in SPM sources

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_MMN')
        
        % Do amplitude thresholds separately for each modality
        for mod = 1:length(modality_data)

            % There will be no average file to obtain MMN from
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
            infolder = find(endsWith({folders.name},[modality_data{mod} '_averaged.mat']) & (...
                        startsWith({folders.name},[participant{p}])));
            if isempty(infolder)
                error(['No ' modality_data{mod} ' averaged files for ' participant{p}]);
            end

            % Should not be by definition, as it is 
            if length(infolder) > 1
                error(['More than one ' modality_data{mod} ' averaged file for ' participant{p}]);
            end

            jobfile = {[subject_file_path '/MMN_job.m']};
            jobs = repmat(jobfile, 1, 1);
            inputs = cell(0, 1);
            inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]);
            eval(['cd ' outcome_path]) % To ensure it saves them there
            spm('defaults', 'EEG');
            spm_jobman('run', jobs, inputs{:});


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust names
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},[modality_data{mod} '_averaged.mat']) & (...
                        startsWith({folders.name},['W' participant{p}])));
            if length(infolder) > 1
                error(['More than one averaged file for ' participant{p} ' ' modality_data{mod}]);
            end      

            % Should only be one, but just in case
            load([outcome_path '/' folders(infolder).name]);
            D.data.fname = [outcome_path '/' participant{p} '_' modality_data{mod} '_MMN_DCM.dat'];
            D.fname = [participant{p} '_' modality_data{mod} '_MMN_DCM.mat'];

            % Adjust name of .dat file with data as well
            datfilename = [outcome_path '/W' participant{p} '_' modality_data{mod} '_averaged.dat'];  %#ok<*SAGROW> 
            % If it was already renamed there will be not imported_dat
            if exist(datfilename,'file')
                datnewname = [outcome_path '/' participant{p} '_' modality_data{mod} '_MMN_DCM.dat'];  %#ok<*SAGROW>  
                movefile(datfilename,datnewname);
            end

            % Now, save .mat (will overwrite previous .mat)
            save([outcome_path '/' participant{p} '_' modality_data{mod} '_MMN_DCM'],'D');   
            % Delete original file
            delete([outcome_path '/W' participant{p} '_' modality_data{mod} '_averaged.mat']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % No need to adjust names (resulting file is named the same)
            disp(' ');      
            disp('-------------------------');
            disp(['MMN DCM file computed for ' participant{p} ' ' modality_data{mod}]);
            disp(datetime)
            disp(' ');                 
            
            % MMN file contains the same than average file plus MMN, so:
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Delete previous step
            if strcmp(delete_average_files,'YES')
                delete([outcome_path '/' participant{p} '_' modality_data{mod} '_averaged.dat']);
                delete([outcome_path '/' participant{p} '_' modality_data{mod} '_averaged.mat']);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        end

        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_forward';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 

    end   
end

%% (Optional) Crop MEGMAG or MEGPLANAR to run in DCM (Only one modality can be run)

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_crop_DCM')
        
        % Do amplitude thresholds separately for each modality
        for mod = 1:length(modality_data)    
            
            % Only if it's MEG
            if strcmp(modality_data{mod},'MEG')
            
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

                jobfile = {[subject_file_path '/Crop_' which_MEG_sensors '_channels_job.m']};
                jobs = repmat(jobfile, 1, 1);
                inputs = cell(0, 1);
                inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]);
                eval(['cd ' outcome_path]) % To ensure it saves them there
                spm('defaults', 'EEG');
                spm_jobman('run', jobs, inputs{:});


                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Adjust names
                folders = dir([outcome_path '/']);
                infolder = find(endsWith({folders.name},'_MMN_DCM.mat') & (...
                            startsWith({folders.name},['C' participant{p} '_' modality_data{mod}])));
                % Overwrite original file name
                for in = 1:length(infolder)
                    line = infolder(in);
                    load([outcome_path '/' folders(line).name]);
                    D.data.fname = [outcome_path '/' participant{p} '_' modality_data{mod} '_MMN_DCM.dat'];
                    D.fname = [participant{p} '_' modality_data{mod} '_MMN_DCM.mat'];

                    % Adjust name of .dat file with data as well
                    datfilename = [outcome_path '/C' participant{p} '_' modality_data{mod} '_MMN_DCM.dat'];  %#ok<*SAGROW> 
                    % If it was already renamed there will be not imported_dat
                    if exist(datfilename,'file')
                        datnewname = [outcome_path '/' participant{p} '_' modality_data{mod} '_MMN_DCM.dat'];  %#ok<*SAGROW>  
                        movefile(datfilename,datnewname);
                    end

                    % Now, save .mat
                    save([outcome_path '/' participant{p} '_' modality_data{mod} '_MMN_DCM'],'D');
                    % Delete original file
                    delete([outcome_path '/C' participant{p} '_' modality_data{mod} '_MMN_DCM.mat']);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % Here never delete the cleaned file

                disp(' ');      
                disp('-------------------------');
                disp([participant{p} ' CROPED DCM ' modality_data{mod}]);
                disp(datetime)
                disp(' '); 
            end 
        end

        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_forward';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 

    end   
end

%% Compute forward model with MMN DCM file

% Will need it even for DCM

% To take into account in SPM forward models (as in manual):
% 1) ~10,000 vertices per hemisphere using the highest (fine) resolution
% 2) The individual cortical mesh is obtained automatically from a canonical mesh in MNI

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_forward')
        
        % Before anything else, ensure MRI exist or copy it in directory
        % If participant folder exists, T1 was probably imported already
        if ~exist([outcome_path '/MRIs/' participant{p} '/' ], 'dir')
            try
                anat_folder_bs = dir(brainstorm_anatomy_path);
                infolder_anat_bs = find(strcmp({anat_folder_bs.name},participant{p}));
                if isempty(infolder_anat_bs)
                    error(['No anatomy folder for ' participant{p}]);
                end
                if length(infolder_anat_bs) > 1
                    error(['More than one anatomy folder for ' participant{p}]);
                end
                % Check which MRI was used in brainstorm
                load([brainstorm_anatomy_path '/' anat_folder_bs(infolder_anat_bs).name '/brainstormsubject.mat'],'Anatomy');
                if isempty(Anatomy) % No Import from string
                    error(['No Anatomy variable for ' participant{p}]);
                end
                load([brainstorm_anatomy_path '/' anat_folder_bs(infolder_anat_bs).name '/' Anatomy(7:end)],'History'); % Because subject names are always 5 characters, plus the '/'
                if isempty(History) % No Import from string
                    error(['No History variable for ' participant{p}]);
                end
                imp_str = find(contains(History,'Import from:'));
                if isempty(imp_str) % No Import from string
                    error(['No import MRI files for ' participant{p}]);
                end
                if length(imp_str) > 1 % more than one Import from string
                    error(['More than 1 import MRI files for ' participant{p}]);
                end
                init_pos = strfind(History{imp_str},'HCPproc');
                T1_name = [HCProc_folder '/' History{imp_str}(init_pos+8:end)];
                [~,~,ext] = fileparts(T1_name);
                if ~strcmp(ext,'.gz')
                    error(['MRI extension cannot be properly decompressed for ' participant{p}]);
                end
                if ~exist(T1_name, 'file')
                    error(['MRI file name in brainstorm is not a file for ' participant{p}]);
                end
                % Create destiny folder to copy it over
                mkdir([outcome_path '/MRIs/', participant{p}]);
                copyfile(T1_name,[outcome_path '/MRIs/' participant{p} '/T1w.nii' ext]);
                
                disp(' ');      
                disp('-------------------------');
                disp(['Extracting .nii for ' participant{p}]);
                disp(datetime)
                disp(' ');    
                gunzip([outcome_path '/MRIs/' participant{p} '/T1w.nii' ext]);
                
            catch
                warning('T1 file does not exist, searching generic one');
                % Try to go for a general one (based on Team's wiki)
                if strcmp(participant{p},'2397A') || strcmp(participant{p},'2229A') % For these two we are using 3mo MRI
                    temp_dir = dir([HCProc_folder '/' participant{p}(1:4) 'B/*/T1w/']);
                    temp_find = find(endsWith({temp_dir.name},'T1w.nii_acpc_dc_restore.nii.gz'));
                    if isempty(temp_find)
                        % Try this instead
                        temp_find = find(endsWith({temp_dir.name},'T1w_acpc_dc_restore.nii.gz'));
                        if isempty(temp_find) % If it is still empty
                            % Try this instead
                            temp_find = find(endsWith({temp_dir.name},'T1w.nii_acpc.nii.gz'));
                            if isempty(temp_find) % If it is still empty
                                error(['nii.acpc.dc.restore.nii.gz file does not exist for ' participant{p}]);
                            end
                        end
                    end
                else
                    temp_dir = dir([HCProc_folder '/' participant{p} '/*/T1w/']);
                    temp_find = find(endsWith({temp_dir.name},'T1w.nii_acpc_dc_restore.nii.gz'));
                    if isempty(temp_find)
                        % Try this instead
                        temp_find = find(endsWith({temp_dir.name},'T1w_acpc_dc_restore.nii.gz'));
                        if isempty(temp_find) % If it is still empty
                            % Try this instead
                            temp_find = find(endsWith({temp_dir.name},'T1w.nii_acpc.nii.gz'));
                            if isempty(temp_find) % If it is still empty
                                error(['nii.acpc.dc.restore.nii.gz file does not exist for ' participant{p}]);
                            end
                        end
                    end
                end
                
                if length(temp_find) > 1 % rare case, but get the final one
                    temp_find = temp_find(end);
                end
                T1_name = [temp_dir(temp_find).folder '/' temp_dir(temp_find).name];  
                [~,~,ext] = fileparts(T1_name);
                if ~strcmp(ext,'.gz')
                    error(['MRI extension cannot be properly decompressed for ' participant{p}]);
                end
                % Create destiny folder to copy it over
                mkdir([outcome_path '/MRIs/', participant{p}]);
                copyfile(T1_name,[outcome_path '/MRIs/' participant{p} '/T1w.nii' ext]);
                
                disp(' ');      
                disp('-------------------------');
                disp(['Extracting .nii for ' participant{p}]);
                disp(datetime)
                disp(' ');    
                gunzip([outcome_path '/MRIs/' participant{p} '/T1w.nii' ext]);
            end
        end
        
        % If we selected fiducials from brainstorm file, calculate them
        if strcmp(fiducials_option,'brainstorm')
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
            mni_bs_coord_nas = cs_convert(sMri, 'mri', 'world', sMri.SCS.NAS ./ 1000) .* 1000;
            mni_bs_coord_lpa = cs_convert(sMri, 'mri', 'world', sMri.SCS.LPA ./ 1000) .* 1000;
            mni_bs_coord_rpa = cs_convert(sMri, 'mri', 'world', sMri.SCS.RPA ./ 1000) .* 1000;
        end
        
        % Compute forward model for all modalities
        for mod = 1:length(modality_data)

            % There will be no MMN file to obtain MMN from
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
            infolder = find(endsWith({folders.name},[modality_data{mod} '_MMN_DCM.mat']) & (...
                        startsWith({folders.name},[participant{p}])));
            if isempty(infolder)
                error(['No ' modality_data{mod} ' MMN DCM files for ' participant{p}]);
            end

            % Should not be by definition, as it is 
            if length(infolder) > 1
                error(['More than one ' modality_data{mod} ' MMN DCM file for ' participant{p}]);
            end

            % If it's a MEG file wih bad EEG channels, forward model on EEG will fail
            if strcmp(modality_data{mod},'MEG') && (strcmp(subject_array{pos_subj,4},'no_EEG_chans') || strcmp(subject_array{pos_subj,5},'bad_EEG'))
               load([outcome_path '/' folders(infolder).name]);
               % Change sensor labels so that it does not take that into account
               for ch = 1:length(D.channels)
                   if strcmp(D.channels(ch).type,'EEG')
                       D.channels(ch).type = 'Other';
                   end
               end
               save([outcome_path '/' folders(infolder).name],'D') 
            end
            
            % If it's an EEG file wih bad MEG channels, forward model on MEG will fail
            if strcmp(modality_data{mod},'EEG') && strcmp(subject_array{pos_subj,6},'bad_MEG')
               load([outcome_path '/' folders(infolder).name]);
               % Change sensor labels so that it does not take that into account
               for ch = 1:length(D.channels)
                   if strcmp(D.channels(ch).type,'MEGMAG') || strcmp(D.channels(ch).type,'MEGPLANAR')
                       D.channels(ch).type = 'Other';
                   end
               end
               save([outcome_path '/' folders(infolder).name],'D') 
            end
                        
            jobfile = {[subject_file_path '/Head_model_job_' fiducials_option '.m']};
            jobs = repmat(jobfile, 1, 1);
            inputs = cell(0, 1);
            % Data file
            inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]); % {'X:\PEPP\User\DCM\SPM_MEG_data\2259A_EEG_MMN.mat'}
            % Forward model files (stored in this same folder)
            inputs{1,2} = {[outcome_path '/MRIs/' participant{p} '/T1w.nii,1']};
            % {'X:\PEPP\User\DCM\SPM_MEG_data\MRIs_uncompressed\2259A_20180605_T1w.nii_acpc_dc_restore.nii,1'}
            % If brainstorm coordinates for fiducials were selected
            if strcmp(fiducials_option,'brainstorm')
                inputs{1,3} = mni_bs_coord_nas;
                inputs{1,4} = mni_bs_coord_lpa;
                inputs{1,5} = mni_bs_coord_rpa;
            end
            eval(['cd ' outcome_path]) % To ensure it saves files there
            spm('defaults', 'EEG');
            spm_jobman('run', jobs, inputs{:});

            % No need to delete any file or change names
            % (as it will include forward in separate file and same file)

            % No need to adjust names (resulting file is named the same)
            disp(' ');      
            disp('-------------------------');
            disp(['Forward model computed for ' participant{p} ' ' modality_data{mod}]);
            disp(datetime)
            disp(' ');                 
            
        end

        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_DCM_CMC';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 

    end   
end
