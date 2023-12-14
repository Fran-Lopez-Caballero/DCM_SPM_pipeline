%% Analyze .fif files EEG/MEG files with SPM (to run DCM)

%% Define variables and start SPM

running_in = 'server'; % 'server' or 'local'

if strcmp(running_in,'server')
    rmpath('/software/path/brainstorm3_v20220706'); % Remove because of overlapping fieltrip functions
    addpath('/software/path/spm12_v7771');
    addpath('/private/path/PEPP/User/DCM/Scripts_final_pipeline');
    subject_file_path = '/private/path/PEPP/User/DCM/Scripts_final_pipeline';
    load('/private/path/PEPP/User/DCM/Scripts_final_pipeline/Labels_EEG.mat');
    root_dir_analysis_path = '/private/path/PEPP/analysis';
    brainstorm_anatomy_path = '/private/path/PEPP/brainstorm_db/PEPP_FEManat_baseline/anat'; % For Baseline
    HCProc_folder = '/private/path/PEPP/HCPproc';
    root_dir_bs = '/private/path/PEPP/brainstorm_db/PEPP_MMN/data';
    outcome_path = '/private/path/PEPP/User/DCM/SPM_MEG_data';
elseif strcmp(running_in,'local')
    rmpath('C:/Users/private/path/private/path/Documents/MATLAB/brainstorm3'); % Remove because of overlapping fieltrip functions
    addpath('C:/Users/private/path/private/path/Documents/MATLAB/spm12');
    addpath('private/path/User/DCM/Scripts_final_pipeline');
    subject_file_path = 'private/path/User/DCM/Scripts_final_pipeline';
    load('private/path/User/DCM/Scripts_final_pipeline/Labels_EEG.mat');
    root_dir_analysis_path = 'private/path/analysis';
    brainstorm_anatomy_path = 'private/path/brainstorm_db/PEPP_FEManat_baseline/anat'; % For Baseline
    HCProc_folder = 'private/path/HCPproc';
    root_dir_bs = 'private/path/brainstorm_db/PEPP_MMN/data';
    outcome_path = 'private/path/User/DCM/SPM_MEG_data';
end

participant_group = {'C','FE'};
subject_array_filename = 'subject_array_matched_baseline';
load([subject_file_path '/' subject_array_filename '.mat']);
modality_data = {'EEG', 'MEG'}; % 'EEG', 'MEG', 'BIMODAL'
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

%% Import AMICA.fif files to SPM

for p = 1:length(participant)

    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_import')
    
    folders = dir([root_dir_analysis_path '/' participant{p} '/']);
    infolder = find(endsWith({folders.name},'tsss_AMICA.fif') & (...
                contains({folders.name},[participant{p} '_MMN'])));
    if isempty(infolder)
        error(['No MMN AMICA files for ' participant{p}]);
    end

    for i = 1:length(infolder)
        line = infolder(i);
        file_name = [root_dir_analysis_path '/' participant{p} '/' folders(line).name];  %#ok<*SAGROW>             

        disp(' ');      
        disp('-------------------------');  
        disp(['Importing EEG/MEG data in SPM for ' participant{p}]);
        disp(datetime)
        disp(' '); 

        S.dataset = file_name;
        S.mode = 'continuous';
        S.outfile = [outcome_path '/' participant{p} '_run' num2str(i) '_imported'];
        S.channels = 'all';
        S.eventpadding = 0;
        S.blocksize = 3276800;
        S.checkboundary = 1;
        S.saveoriginheader = 0;
        [D] = spm_eeg_convert(S); % Will save files in local directory
    end
    
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_scale';
    save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
    
    disp(' ');      
    disp('-------------------------');
    disp([participant{p} ' IMPORTED']);
    disp(datetime)
    disp(' '); 
    
    end
end

%% Scale MEG if needed

% Because of some subjects having wrong MEG scales (preprocessing problem
% in 2021...) we need to adjust that by applying a montage
% Since SPM will always leave a prefix in files and to keep consistency, we
% will scale all files, only that the ones that don't need it will be
% multiplied by 1 to keep their values the same

for p = 1:length(participant)

    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p})); %#ok<*CCAT1>
    if strcmp(subject_array{pos_subj,3},'needs_scale') && strcmp(subject_array{pos_subj,4},'scale')
    
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_imported.mat') & (...
                    startsWith({folders.name},[participant{p}])));
        if isempty(infolder)
            error(['No imported files for ' participant{p}]);
        end

        % List of open inputs
        nrun = length(infolder); % enter the number of runs here
        jobfile = {[subject_file_path '/Scale_job_' running_in '.m']};
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(0, nrun);
        for crun = 1:nrun
            line = infolder(crun);
            inputs{1,crun} = cellstr([outcome_path '/' folders(line).name]);
        end
        spm('defaults', 'EEG');
        spm_jobman('run', jobs, inputs{:});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adjust names
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_imported.mat') & (...
                    startsWith({folders.name},['S' participant{p}])));
        % Rename each run file
        for in = 1:length(infolder)
            line = infolder(in);
            load([outcome_path '/' folders(line).name]);
            D.data.fname = [outcome_path '/' participant{p} '_run' num2str(in) '_scaled.dat'];
            D.fname = [participant{p} '_run' num2str(in) '_scaled.mat'];
            
            % Adjust name of .dat file with data as well
            datfilename = [outcome_path '/S' participant{p} '_run' num2str(in) '_imported.dat'];  %#ok<*SAGROW> 
            % If it was already renamed there will be not imported_dat
            if exist(datfilename,'file')
                datnewname = [outcome_path '/' participant{p} '_run' num2str(in) '_scaled.dat'];  %#ok<*SAGROW>  
                movefile(datfilename,datnewname);
            end
            % Save
            save([outcome_path '/' participant{p} '_run' num2str(in) '_scaled'],'D');
            % Delete original file
            delete([outcome_path '/S' participant{p} '_run' num2str(in) '_imported.mat']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_correction';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Delete previous step
        if strcmp(delete_previous_steps,'YES')
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_imported.mat') & (...
                        startsWith({folders.name},[participant{p}])));
            for in = 1:length(infolder)
                % Adjust name of .dat file with data as well
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_imported.dat']);
                % Delete original file
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_imported.mat']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp(' ');      
        disp('-------------------------');
        disp([participant{p} ' SCALED (MEG)']);
        disp(datetime)
        disp(' '); 
        
    elseif strcmp(subject_array{pos_subj,3},'needs_scale') && ~strcmp(subject_array{pos_subj,4},'scale')
    
        % Adjust names
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_imported.mat') & (...
                    startsWith({folders.name},[participant{p}])));
        % Rename each run file
        for in = 1:length(infolder)
            line = infolder(in);
            load([outcome_path '/' folders(line).name]);
            D.data.fname = [outcome_path '/' participant{p} '_run' num2str(in) '_scaled.dat'];
            D.fname = [participant{p} '_run' num2str(in) '_scaled.mat'];
            
            % Adjust name of .dat file with data as well
            datfilename = [outcome_path '/' participant{p} '_run' num2str(in) '_imported.dat'];  %#ok<*SAGROW> 
            % If it was already renamed there will be not imported_dat
            if exist(datfilename,'file')
                datnewname = [outcome_path '/' participant{p} '_run' num2str(in) '_scaled.dat'];  %#ok<*SAGROW>  
                copyfile(datfilename,datnewname);
            end
            % Save
            save([outcome_path '/' participant{p} '_run' num2str(in) '_scaled'],'D');
            % Delete original file
            % delete([outcome_path '/' participant{p} '_run' num2str(in) '_imported.mat']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_correction';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Delete previous step
        if strcmp(delete_previous_steps,'YES')
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_imported.mat') & (...
                        startsWith({folders.name},[participant{p}])));
            for in = 1:length(infolder)
                % Adjust name of .dat file with data as well
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_imported.dat']);
                % Delete original file
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_imported.mat']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        disp(' ');      
        disp('-------------------------');
        disp([participant{p} ' SCALED (MEG)']);
        disp(datetime)
        disp(' '); 
    end
end

%% Correct epochs and channels based on brainstorm file

for p = 1:length(participant)

    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_correction')
    
    folders = dir([outcome_path '/']);
    infolder = find(endsWith({folders.name},'_scaled.mat') & (...
                    startsWith({folders.name},[participant{p}])));
    
    if isempty(infolder)
        error(['No SPM scaled files for ' participant{p}]);
    end

    for i = 1:length(infolder) % for each run if there is more than one
        line = infolder(i);
        file_name = [outcome_path '/' folders(line).name];  %#ok<*SAGROW>             

        disp(' ');      
        disp('-------------------------');  
        disp(['Correcting events and channels for ' participant{p}]);
        disp(datetime)
        disp(' '); 

        % Load 'D' to be modified
        load(file_name);

        % Modify triggers
        if strcmp(modify_triggers,'YES')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load vmrk file of the corresponding run
        delimiter = ',';
        startRow = 12;
        formatSpec = '%*q%q%q%*s%*s%*s%[^\n\r]';
        foldersvmrk = dir([root_dir_analysis_path '/' participant{p} '/']);
        if length(infolder) > 1 % there is more than one run
            infoldervmrk = find(endsWith({foldersvmrk.name},'tsss_AMICA.vmrk') & (...
            contains({foldersvmrk.name},[participant{p} '_MMN_run' num2str(i)])));
        else
            infoldervmrk = find(endsWith({foldersvmrk.name},'tsss_AMICA.vmrk') & (...
            contains({foldersvmrk.name},[participant{p} '_MMN_'])));
        end

        if isempty(infoldervmrk)
            error(['No vmrk files for ' participant{p} '_run' num2str(i)]);
        end
        if length(infoldervmrk) > 1
            error(['More than one vmrk files for ' participant{p} '_run' num2str(i)]);
        end
        % linevmrk = infoldervmrk(i);
        fullfilename = [root_dir_analysis_path '/' participant{p} '/' foldersvmrk(infoldervmrk).name];  %#ok<*SAGROW>   
        fileID = fopen(fullfilename,'r');
        textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'ReturnOnError', false);
        fclose(fileID);
        
        % Update 'D.trials.events' based on 
        replacement_triggers = D.trials.events(1:2); % For structure
        % Replace strings with labels
        for l = 1:length(dataArray{1,1})
            if ~isempty(dataArray{1,1}{l})
                replacement_triggers(l).type = dataArray{1,1}{l};
            else
                replacement_triggers(l).type = 'empty';
            end
        end
        % Replace values (in case this indicates triggers
        for l = 1:length(dataArray{1,1})
            if ~isnan(str2double(dataArray{1,1}{l}))
                replacement_triggers(l).value = str2double(dataArray{1,1}{l});
            else % It's a boundary or empty
                replacement_triggers(l).value = 0;
            end
        end
        % Replace times
        for l = 1:length(dataArray{1,2})
            % 8032 to 8.0320
            replacement_triggers(l).time = (str2double(dataArray{1,2}{l}))/1000;
        end
        % Correct offset
        for l = 1:length(dataArray{1,1})
            replacement_triggers(l).offset = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
        % Modify EEG channels
        if ~strcmp(subject_array{pos_subj,4},'no_EEG_chans')
            % Now, correct EEG channel labels (MEG ones are ok)
            % D.sensors.eeg.label
            index_to_keep = []; pos = 1;
            for ei = 1:length(Equiv_table)
                index_elec = find(strcmp(D.sensors.eeg.label,Equiv_table{ei,1}));
                if ~isempty(index_elec)
                    D.sensors.eeg.label{index_elec} = Equiv_table{ei,2};
                    index_to_keep(pos) = index_elec;
                    pos = pos + 1;
                else
                    error(['Electrode ' Equiv_table{ei,1} 'not found for ' participant{p}]);
                end
            end

            % Now remove any extra electrodes      

            % Same for D.channels:
            index_to_keep = []; pos = 1;
            for ei = 1:length(Equiv_table)
                index_elec = find(strcmp({D.channels.label},Equiv_table{ei,1}));
                if ~isempty(index_elec)
                    D.channels(index_elec).label = Equiv_table{ei,2};
                    % Also, ensure this says EEG:
                    D.channels(index_elec).type = 'EEG';
                    index_to_keep(pos) = index_elec;
                    pos = pos + 1;
                else
                    error(['Electrode ' Equiv_table{ei,1} 'not found in channel.label for ' participant{p}]);
                end
            end

            % Manually replace EEG type of EEG 61 to 64
            try
                D.channels(find(strcmp({D.channels.label},'EEG061'))).type = 'Other'; %#ok<*FNDSB>
            catch
                D.channels(find(strcmp({D.channels.label},'EOG061'))).type = 'Other'; %#ok<*FNDSB>
            end
            try 
                D.channels(find(strcmp({D.channels.label},'EEG062'))).type = 'Other'; %#ok<*FNDSB>
            catch
                D.channels(find(strcmp({D.channels.label},'EOG062'))).type = 'Other'; %#ok<*FNDSB>
            end
            try
                D.channels(find(strcmp({D.channels.label},'EEG063'))).type = 'Other'; %#ok<*FNDSB>
            catch
                D.channels(find(strcmp({D.channels.label},'EOG063'))).type = 'Other'; %#ok<*FNDSB>
            end
            try
                D.channels(find(strcmp({D.channels.label},'EEG064'))).type = 'Other'; %#ok<*FNDSB>
            catch
                D.channels(find(strcmp({D.channels.label},'EOG064'))).type = 'Other'; %#ok<*FNDSB>
            end
        end
        
        % Add bad channels from brainstorm file
        % Good bad channels are the same across conditions, so Standard
        foldersbs = dir([root_dir_bs '/' participant{p} '/Standard']);
        infolderbs = find(contains({foldersbs.name},'trial'));
        if isempty(infolderbs)
            error(['no brainstorm standard trial data for ' participant{p}]);
        else
            line = infolderbs(1); % any trial will do
            load([root_dir_bs '/' participant{p} '/Standard/' foldersbs(line).name]);
        end
        try
            load([root_dir_bs '/' participant{p} '/Standard/channel_vectorview306_acc1.mat']);
        catch
            error(['no brainstorm standard channel file for ' participant{p}]);
        end

        for cf = 1:length(D.channels)
            if ~strcmp(D.channels(cf).type,'EEG') && ~strcmp(D.channels(cf).type,'MEGPLANAR') && ~strcmp(D.channels(cf).type,'MEGMAG')
                continue
            end

            if strcmp(subject_array{pos_subj,4},'no_EEG_chans') && strcmp(D.channels(cf).type,'EEG')
                continue
            end
            pos_bc = find(strcmp({Channel.Name},D.channels(cf).label));
            if isempty(pos_bc)
                error(['no ' D.channels(cf).label 'found in bs for ' participant{p}]);
            else
                if ChannelFlag(pos_bc) == 1
                    D.channels(cf).bad = 0;
                elseif ChannelFlag(pos_bc) == -1
                    D.channels(cf).bad = 1;
                end
            end
        end

        % Adjust names so that it can be read by SPM
        D.data.fname = [outcome_path '/' participant{p} '_run' num2str(i) '_corrected.dat'];
        D.fname = [participant{p} '_run' num2str(i) '_corrected.mat'];
        % Adjust name of .dat file with data as well

        datfilename = [outcome_path '/' participant{p} '_run' num2str(i) '_scaled.dat'];  %#ok<*SAGROW> 
        % If it was already renamed there will be not imported_dat
        if exist(datfilename,'file')
            datnewname= [outcome_path '/' participant{p} '_run' num2str(i) '_corrected.dat'];  %#ok<*SAGROW>  
            copyfile(datfilename,datnewname);
        end
        % Save
        if strcmp(modify_triggers,'YES')
            D.trials.events = replacement_triggers;
        end
        save([outcome_path '/' participant{p} '_run' num2str(i) '_corrected'],'D');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Delete previous step
    if strcmp(delete_previous_steps,'YES')
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_scaled.mat') & (...
                    startsWith({folders.name},[participant{p}])));
        for in = 1:length(infolder)
            % Adjust name of .dat file with data as well
            delete([outcome_path '/' participant{p} '_run' num2str(in) '_scaled.dat']);
            % Delete original file
            delete([outcome_path '/' participant{p} '_run' num2str(in) '_scaled.mat']);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % If successful, update subject_array for this subject
    subject_array{pos_subj,3} = 'needs_reref';
    save([subject_file_path '/' subject_array_filename '.mat'],'subject_array')
    
    disp(' ');      
    disp('-------------------------');
    disp([participant{p} ' CORRECTED']);
    disp(datetime)
    disp(' '); 
    end
end

%% Re-reference EEG

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_reref') && ~strcmp(subject_array{pos_subj,4},'no_EEG_chans') 
        % ... && ~strcmp(subject_array{pos_subj,5},'bad_EEG') % Has EEG channels to do reref
    
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_corrected.mat') & (...
                    startsWith({folders.name},[participant{p}])));
        if isempty(infolder)
            error(['No corrected files for ' participant{p}]);
        end

        % List of open inputs
        nrun = length(infolder); % enter the number of runs here
        jobfile = {[subject_file_path '/EEG_avg_reref_job.m']};
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(0, nrun);
        for crun = 1:nrun
            line = infolder(crun);
            inputs{1,crun} = cellstr([outcome_path '/' folders(line).name]);
        end
        spm('defaults', 'EEG');
        spm_jobman('run', jobs, inputs{:});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adjust names
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_corrected.mat') & (...
                    startsWith({folders.name},['M' participant{p}])));
        % Rename each run file
        for in = 1:length(infolder)
            line = infolder(in);
            load([outcome_path '/' folders(line).name]);
            D.data.fname = [outcome_path '/' participant{p} '_run' num2str(in) '_reref.dat'];
            D.fname = [participant{p} '_run' num2str(in) '_reref.mat'];
            
            % Adjust name of .dat file with data as well
            datfilename = [outcome_path '/M' participant{p} '_run' num2str(in) '_corrected.dat'];  %#ok<*SAGROW> 
            % If it was already renamed there will be not imported_dat
            if exist(datfilename,'file')
                datnewname = [outcome_path '/' participant{p} '_run' num2str(in) '_reref.dat'];  %#ok<*SAGROW>  
                movefile(datfilename,datnewname);
            end
            
            % When re-referencing EEG, for some reason EEG channels are set to
            % 'Other' type and lose positions in channel variable (not in
            % sensor one): so we correct this by copying the one from the
            % _corrected.mat file
            
            % Temporarily save original D
            D_orig = D;
            % load D from _corrected.mat
            load([outcome_path '/' participant{p} '_run' num2str(in) '_corrected.mat']);
            % D_orig.channels
            pos_eeg_corr = find(strcmp({D.channels.type},'EEG'));
            eeg_chan_labels = {D.channels(pos_eeg_corr).label};
            for pec = 1:length(eeg_chan_labels)
                pos_cc = find(strcmp({D.channels.label},eeg_chan_labels{pec}));
                pos_co = find(strcmp({D_orig.channels.label},eeg_chan_labels{pec}));
                D_orig.channels(pos_co).type = D.channels(pos_cc).type;
                D_orig.channels(pos_co).X_plot2D = D.channels(pos_cc).X_plot2D;
                D_orig.channels(pos_co).Y_plot2D = D.channels(pos_cc).Y_plot2D;
                % These next two don't seem to change but just in case
                D_orig.channels(pos_co).bad = D.channels(pos_cc).bad;
                D_orig.channels(pos_co).units = D.channels(pos_cc).units;
            end
            % Overwrite D with the D_orig that was corrected
            D = D_orig;
            % Now, save .mat
            save([outcome_path '/' participant{p} '_run' num2str(in) '_reref'],'D');
            % Delete original file
            delete([outcome_path '/M' participant{p} '_run' num2str(in) '_corrected.mat']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Delete previous step
        if strcmp(delete_previous_steps,'YES')
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_corrected.mat') & (...
                        startsWith({folders.name},[participant{p}])));
            for in = 1:length(infolder)
                % Adjust name of .dat file with data as well
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_corrected.dat']);
                % Delete original file
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_corrected.mat']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_filter';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
        
        disp(' ');      
        disp('-------------------------');
        disp([participant{p} ' REREFERED (EEG)']);
        disp(datetime)
        disp(' '); 
    
    elseif strcmp(subject_array{pos_subj,3},'needs_reref') && strcmp(subject_array{pos_subj,4},'no_EEG_chans') % Does not have EEG channels
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adjust names
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_corrected.mat') & (...
                    startsWith({folders.name},[participant{p}])));
        % Rename each run file
        for in = 1:length(infolder)
            line = infolder(in);
            load([outcome_path '/' folders(line).name]);
            D.data.fname = [outcome_path '/' participant{p} '_run' num2str(in) '_reref.dat'];
            D.fname = [participant{p} '_run' num2str(in) '_reref.mat'];
            
            % Adjust name of .dat file with data as well
            datfilename = [outcome_path '/' participant{p} '_run' num2str(in) '_corrected.dat'];  %#ok<*SAGROW> 
            % If it was already renamed there will be not imported_dat
            if exist(datfilename,'file')
                datnewname = [outcome_path '/' participant{p} '_run' num2str(in) '_reref.dat'];  %#ok<*SAGROW>  
                copyfile(datfilename,datnewname);
            end
            % Save
            save([outcome_path '/' participant{p} '_run' num2str(in) '_reref'],'D');
            % Delete original file
            % delete([outcome_path '/M' participant{p} '_run' num2str(in) '_corrected.mat']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Delete previous step
        if strcmp(delete_previous_steps,'YES')
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_corrected.mat') & (...
                        startsWith({folders.name},[participant{p}])));
            for in = 1:length(infolder)
                % Adjust name of .dat file with data as well
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_corrected.dat']);
                % Delete original file
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_corrected.mat']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_filter';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array')
    
        disp(' ');      
        disp('-------------------------');
        disp([participant{p} ' REREFERED (EEG)']);
        disp(datetime)
        disp(' '); 
    end   
end

%% Filter

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_filter')
    
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_reref.mat') & (...
                    startsWith({folders.name},[participant{p}])));
        if isempty(infolder)
            error(['No reref files for ' participant{p}]);
        end

        % List of open inputs
        nrun = length(infolder); % enter the number of runs here
        jobfile = {[subject_file_path '/Filter_job.m']};
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(0, nrun);
        for crun = 1:nrun
            line = infolder(crun);
            inputs{1,crun} = cellstr([outcome_path '/' folders(line).name]);
        end
        spm('defaults', 'EEG');
        spm_jobman('run', jobs, inputs{:});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adjust names
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_reref.mat') & (...
                    startsWith({folders.name},['F' participant{p}])));
        % Rename each run file
        for in = 1:length(infolder)
            line = infolder(in);
            load([outcome_path '/' folders(line).name]);
            D.data.fname = [outcome_path '/' participant{p} '_run' num2str(in) '_filtered.dat'];
            D.fname = [participant{p} '_run' num2str(in) '_filtered.mat'];
            
            % Adjust name of .dat file with data as well
            datfilename = [outcome_path '/F' participant{p} '_run' num2str(in) '_reref.dat'];  %#ok<*SAGROW> 
            % If it was already renamed there will be not imported_dat
            if exist(datfilename,'file')
                datnewname = [outcome_path '/' participant{p} '_run' num2str(in) '_filtered.dat'];  %#ok<*SAGROW>  
                movefile(datfilename,datnewname);
            end
            
            % Now, save .mat
            save([outcome_path '/' participant{p} '_run' num2str(in) '_filtered'],'D');
            % Delete original file
            delete([outcome_path '/F' participant{p} '_run' num2str(in) '_reref.mat']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Delete previous step
        if strcmp(delete_previous_steps,'YES')
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_reref.mat') & (...
                        startsWith({folders.name},[participant{p}])));
            for in = 1:length(infolder)
                % Adjust name of .dat file with data as well
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_reref.dat']);
                % Delete original file
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_reref.mat']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_epoch';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
    
        disp(' ');      
        disp('-------------------------');
        disp([participant{p} ' FILTERED']);
        disp(datetime)
        disp(' '); 
        
    end   
end

%% Epoch and baseline correction

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_epoch')
    
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_filtered.mat') & (...
                    startsWith({folders.name},[participant{p}])));
        if isempty(infolder)
            error(['No reref files for ' participant{p}]);
        end

        % List of open inputs
        nrun = length(infolder); % enter the number of runs here
        if strcmp(modify_triggers,'YES')
            jobfile = {[subject_file_path '/Epoch_new_triggers_job.m']};
        elseif strcmp(modify_triggers,'NO')
            jobfile = {[subject_file_path '/Epoch_original_triggers_job.m']};
        end
        jobs = repmat(jobfile, 1, nrun);
        inputs = cell(0, nrun);
        for crun = 1:nrun
            line = infolder(crun);
            inputs{1,crun} = cellstr([outcome_path '/' folders(line).name]);
        end
        spm('defaults', 'EEG');
        spm_jobman('run', jobs, inputs{:});
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Adjust names
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_filtered.mat') & (...
                    startsWith({folders.name},['E' participant{p}])));
        % Rename each run file
        for in = 1:length(infolder)
            line = infolder(in);
            load([outcome_path '/' folders(line).name]);
            D.data.fname = [outcome_path '/' participant{p} '_run' num2str(in) '_epoched.dat'];
            D.fname = [participant{p} '_run' num2str(in) '_epoched.mat'];
            
            % Adjust name of .dat file with data as well
            datfilename = [outcome_path '/E' participant{p} '_run' num2str(in) '_filtered.dat'];  %#ok<*SAGROW> 
            % If it was already renamed there will be not imported_dat
            if exist(datfilename,'file')
                datnewname = [outcome_path '/' participant{p} '_run' num2str(in) '_epoched.dat'];  %#ok<*SAGROW>  
                movefile(datfilename,datnewname);
            end
            
            % Now, save .mat
            save([outcome_path '/' participant{p} '_run' num2str(in) '_epoched'],'D');
            % Delete original file
            delete([outcome_path '/E' participant{p} '_run' num2str(in) '_filtered.mat']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Delete previous step
        if strcmp(delete_previous_steps,'YES')
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_filtered.mat') & (...
                        startsWith({folders.name},[participant{p}])));
            for in = 1:length(infolder)
                % Adjust name of .dat file with data as well
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_filtered.dat']);
                % Delete original file
                delete([outcome_path '/' participant{p} '_run' num2str(in) '_filtered.mat']);
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_merge';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
        
        disp(' ');      
        disp('-------------------------');
        disp([participant{p} ' EPOCHED']);
        disp(datetime)
        disp(' '); 
    
    end   
end

%% Merge runs

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_merge')
    
        folders = dir([outcome_path '/']);
        infolder = find(endsWith({folders.name},'_epoched.mat') & (...
                    startsWith({folders.name},[participant{p}])));
        if isempty(infolder)
            error(['No reref files for ' participant{p}]);
        end

        if length(infolder) > 1
            % Create cell array with names to feed the SPM job
            files_to_merge = {};
            for j = 1:length(infolder)
                whichline = infolder(j);
                files_to_merge{j,1} = [outcome_path '/' folders(whichline).name];
            end
            
            % List of open inputs
            jobfile = {[subject_file_path '/Merge_job.m']};
            jobs = repmat(jobfile, 1, 1);
            inputs = cell(0, 1);
            inputs{1,1} = files_to_merge;
            % Change present path as, for some reason, it saves it  there
            eval(['cd ' outcome_path])
            spm('defaults', 'EEG');
            spm_jobman('run', jobs, inputs{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust names
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_epoched.mat') & (...
                        startsWith({folders.name},['C' participant{p}])));
            % Should only be one, but just in case
            for in = 1:length(infolder)
                line = infolder(in);
                load([outcome_path '/' folders(line).name]);
                D.data.fname = [outcome_path '/' participant{p} '_merged.dat'];
                D.fname = [participant{p} '_merged.mat'];

                % Adjust name of .dat file with data as well
                datfilename = [outcome_path '/C' participant{p} '_run' num2str(in) '_epoched.dat'];  %#ok<*SAGROW> 
                % If it was already renamed there will be not imported_dat
                if exist(datfilename,'file')
                    datnewname = [outcome_path '/' participant{p} '_merged.dat'];  %#ok<*SAGROW>  
                    movefile(datfilename,datnewname);
                end

                % Now, save .mat
                save([outcome_path '/' participant{p} '_merged'],'D');
                % Delete original file
                delete([outcome_path '/C' participant{p} '_run' num2str(in) '_epoched.mat']);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Delete previous step
            if strcmp(delete_previous_steps,'YES')
                folders = dir([outcome_path '/']);
                infolder = find(endsWith({folders.name},'_epoched.mat') & (...
                            startsWith({folders.name},[participant{p}])));
                for in = 1:length(infolder)
                    % Adjust name of .dat file with data as well
                    delete([outcome_path '/' participant{p} '_run' num2str(in) '_epoched.dat']);
                    % Delete original file
                    delete([outcome_path '/' participant{p} '_run' num2str(in) '_epoched.mat']);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % If successful, update subject_array for this subject
            subject_array{pos_subj,3} = 'needs_cleaning';
            save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
            
            disp(' ');      
            disp('-------------------------');
            disp([participant{p} ' MERGED']);
            disp(datetime)
            disp(' '); 
        else
            % Do not merge but change names
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust names
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_epoched.mat') & (...
                        startsWith({folders.name},[participant{p}])));
            % By definition this should not happen, but just in case
            if length(infolder) > 1
                error(['More than one run but not merging for ' participant{p}]);
            end
            
            % Should only be one, but just in case
            for in = 1:length(infolder)
                line = infolder(in);
                load([outcome_path '/' folders(line).name]);
                D.data.fname = [outcome_path '/' participant{p} '_merged.dat'];
                D.fname = [participant{p} '_merged.mat'];

                % Adjust name of .dat file with data as well
                datfilename = [outcome_path '/' participant{p} '_run' num2str(in) '_epoched.dat'];  %#ok<*SAGROW> 
                % If it was already renamed there will be not imported_dat
                if exist(datfilename,'file')
                    datnewname = [outcome_path '/' participant{p} '_merged.dat'];  %#ok<*SAGROW>  
                    copyfile(datfilename,datnewname);
                end

                % Now, save .mat
                save([outcome_path '/' participant{p} '_merged'],'D');
                % Delete original file
                % delete([outcome_path '/C' participant{p} '_run' num2str(in) '_epoched.mat']);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Delete previous step
            if strcmp(delete_previous_steps,'YES')
                folders = dir([outcome_path '/']);
                infolder = find(endsWith({folders.name},'_epoched.mat') & (...
                            startsWith({folders.name},[participant{p}])));
                for in = 1:length(infolder)
                    % Adjust name of .dat file with data as well
                    delete([outcome_path '/' participant{p} '_run' num2str(in) '_epoched.dat']);
                    % Delete original file
                    delete([outcome_path '/' participant{p} '_run' num2str(in) '_epoched.mat']);
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % If successful, update subject_array for this subject
            subject_array{pos_subj,3} = 'needs_cleaning';
            save([subject_file_path '/' subject_array_filename '.mat'],'subject_array')
            
            disp(' ');      
            disp('-------------------------');
            disp([participant{p} ' MERGED']);
            disp(datetime)
            disp(' '); 
        end
    end   
end

%% Amplitude threshold

% * Bear in mind that events too close to the boundaries (and thus containing
% a bunch of "zero" values) could be in the data, so cleaning job in this
% section already includes a flat segment threshold (epochs with more than
% 4 samples with the same value (most likely zero) are discarded

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_cleaning')

        % Define cell arrays to store bad trial count
        sweep_numbers = {};
        sweep_numbers{1,1} = 'Conditions';
        sweep_numbers{1,2} = [participant{p} '_EEG'];
        sweep_numbers{1,3} = [participant{p} '_MEG'];
        sweep_numbers{1,4} = [participant{p} '_BIMODAL'];
        sweep_numbers{2,1} = 'STD_ori';
        sweep_numbers{3,1} = 'STD_thr';
        sweep_numbers{4,1} = 'STD_%';
        sweep_numbers{5,1} = 'PDev_ori';
        sweep_numbers{6,1} = 'PDev_thr';
        sweep_numbers{7,1} = 'PDev_%';
        sweep_numbers{8,1} = 'DDev_ori';
        sweep_numbers{9,1} = 'DDev_thr';
        sweep_numbers{10,1} = 'DDev_%';
        
        % Do amplitude thresholds separately for each modality
        for mod = 1:length(modality_data)
            
            % If there are no EEG channels, don't apply this
            if strcmp(modality_data{mod},'EEG') && strcmp(subject_array{pos_subj,4},'no_EEG_chans')
                continue;
            end
            
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_merged.mat') & (...
                        startsWith({folders.name},[participant{p}])));
            if isempty(infolder)
                error(['No merged files for ' participant{p}]);
            end

            % Should not be by definition, as it is 
            if length(infolder) > 1
                error(['More than one merged file for ' participant{p}]);
            end
            
            jobfile = {[subject_file_path '/' modality_data{mod} '_cleaning_job.m']};
            jobs = repmat(jobfile, 1, 1);
            inputs = cell(0, 1);
            inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]);
            eval(['cd ' outcome_path]) % To ensure it saves them there
            spm('defaults', 'EEG');
            spm_jobman('run', jobs, inputs{:});

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Adjust names
            folders = dir([outcome_path '/']);
            infolder = find(endsWith({folders.name},'_merged.mat') & (...
                        startsWith({folders.name},['A' participant{p}])));
            if length(infolder) > 1
                error(['More than one cleaned file for ' participant{p}]);
            end      
            
            % Should only be one, but just in case
            load([outcome_path '/' folders(infolder).name]);
            D.data.fname = [outcome_path '/' participant{p} '_' modality_data{mod} '_cleaned.dat'];
            D.fname = [participant{p} '_' modality_data{mod} '_cleaned.mat'];

            % Adjust name of .dat file with data as well
            datfilename = [outcome_path '/A' participant{p} '_merged.dat'];  %#ok<*SAGROW> 
            % If it was already renamed there will be not imported_dat
            if exist(datfilename,'file')
                datnewname = [outcome_path '/' participant{p} '_' modality_data{mod} '_cleaned.dat'];  %#ok<*SAGROW>  
                movefile(datfilename,datnewname);
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find total number of sweeps
            trials_count_std = 0;
            trials_count_pdev = 0;
            trials_count_ddev = 0;
            for t = 1:length(D.trials)
                if strcmp(D.trials(t).label, 'Standard')
                    trials_count_std = trials_count_std + 1;
                elseif strcmp(D.trials(t).label, 'PDev')
                    trials_count_pdev = trials_count_pdev + 1;
                elseif strcmp(D.trials(t).label, 'DDev')
                    trials_count_ddev = trials_count_ddev + 1;
                end
            end
            % Find number of surivivng sweeps
            bad_trials_count_std = 0;
            bad_trials_count_pdev = 0;
            bad_trials_count_ddev = 0;
            for t = 1:length(D.trials)
                if D.trials(t).bad == 1 && strcmp(D.trials(t).label, 'Standard')
                    bad_trials_count_std = bad_trials_count_std + 1;
                elseif D.trials(t).bad == 1 && strcmp(D.trials(t).label, 'PDev')
                    bad_trials_count_pdev = bad_trials_count_pdev + 1;
                elseif D.trials(t).bad == 1 && strcmp(D.trials(t).label, 'DDev')
                    bad_trials_count_ddev = bad_trials_count_ddev + 1;
                end
            end
            percentage_std = ((trials_count_std-bad_trials_count_std)/trials_count_std)*100;
            percentage_pdev = ((trials_count_pdev-bad_trials_count_pdev)/trials_count_pdev)*100;
            percentage_ddev = ((trials_count_ddev-bad_trials_count_ddev)/trials_count_ddev)*100;
            surviving_std = trials_count_std-bad_trials_count_std;
            surviving_pdev = trials_count_pdev-bad_trials_count_pdev;
            surviving_ddev = trials_count_ddev-bad_trials_count_ddev;
            
            sweep_numbers{2,mod+1} = num2str(trials_count_std); % 'STD_ori'
            sweep_numbers{3,mod+1} = num2str(surviving_std); % 'STD_thr'
            sweep_numbers{4,mod+1} = [num2str(percentage_std) '%']; % 'STD_%'
            sweep_numbers{5,mod+1} = num2str(trials_count_pdev); % 'PDev_ori'
            sweep_numbers{6,mod+1} = num2str(surviving_pdev); % 'PDev_thr'
            sweep_numbers{7,mod+1} = [num2str(percentage_pdev) '%']; % 'PDev_%'
            sweep_numbers{8,mod+1} = num2str(trials_count_ddev); % 'DDev_ori'
            sweep_numbers{9,mod+1} = num2str(surviving_ddev); % 'DDev_thr'
            sweep_numbers{10,mod+1} = [num2str(percentage_ddev) '%']; % 'DDev_%'
            
            % Take the chance to determine data quality
            if (percentage_std < crit_percent) || (surviving_std < crit_sweeps) || (percentage_pdev < crit_percent) || (surviving_pdev < crit_sweeps) || (percentage_ddev < crit_percent) || (surviving_ddev < crit_sweeps)
                if strcmp(modality_data{mod},'EEG')
                    if ~strcmp(subject_array{pos_subj,5},'exception_EEG')
                        subject_array{pos_subj,5} = 'bad_EEG';
                    end
                elseif strcmp(modality_data{mod},'MEG')
                    if ~strcmp(subject_array{pos_subj,6},'exception_MEG')
                        subject_array{pos_subj,6} = 'bad_MEG';
                    end
                elseif strcmp(modality_data{mod},'BIMODAL')
                    if ~strcmp(subject_array{pos_subj,7},'exception_BIMODAL')
                        subject_array{pos_subj,7} = 'bad_BIMODAL';
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Now, save .mat
            save([outcome_path '/' participant{p} '_' modality_data{mod} '_cleaned'],'D');
            % Delete original file
            delete([outcome_path '/A' participant{p} '_merged.mat']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            disp(' ');      
            disp('-------------------------');
            disp([participant{p} ' CLEANED ' modality_data{mod}]);
            disp(datetime)
            disp(' '); 
        end
        
        % Save sweep count
        save([outcome_path '/QC/' participant{p} '_sweeps'],'sweep_numbers');
        
        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_sorting';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Delete previous step
        if strcmp(delete_merged_files,'YES')
            delete([outcome_path '/' participant{p} '_merged.dat']);
            delete([outcome_path '/' participant{p} '_merged.mat']);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end   
end

%% Sort conditions

% !It is very important so that MMN operations are correct across subjects
% * Refer to manual for specificts, but it's just another "housekeeping"
% cleaning step: file:///C:/Users/private/path/Downloads/fnins-13-00300.pdf

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_sorting')
        
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
            
            jobfile = {[subject_file_path '/Sort_job.m']};
            jobs = repmat(jobfile, 1, 1);
            inputs = cell(0, 1);
            inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]);
            eval(['cd ' outcome_path]) % To ensure it saves them there
            spm('defaults', 'EEG');
            spm_jobman('run', jobs, inputs{:});

            % No need to adjust names (resulting file is named the same)

            disp(' ');      
            disp('-------------------------');
            disp([participant{p} ' SORTED ' modality_data{mod}]);
            disp(datetime)
            disp(' '); 
        end

        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_planar';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
        
        % Do not delete previous steps as the resulting file is the same

    end   
end

%% (OPTIONAL) Combine planar gradiometers (for MEG sensor only)

% It creates "extra" combined MEG sensors that make further steps be longer

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_planar')
        
        % Do amplitude thresholds separately for each modality
        for mod = 1:length(modality_data)
            
            % Only if it is MEG or BIMODAL
            if (strcmp(modality_data{mod},'MEG') && strcmp(combine_planar,'YES') && ~strcmp(subject_array{pos_subj,6},'bad_MEG')) ...
                    || (strcmp(modality_data{mod},'BIMODAL') && strcmp(combine_planar,'YES') && ~strcmp(subject_array{pos_subj,7},'bad_BIMODAL'))
                
                folders = dir([outcome_path '/']);
                infolder = find(endsWith({folders.name},[modality_data{mod} '_cleaned.mat']) & (...
                            startsWith({folders.name},[participant{p}])));
                if isempty(infolder)
                    error(['No ' modality_data{mod} ' merged files for ' participant{p}]);
                end

                % Should not be by definition, as it is 
                if length(infolder) > 1
                    error(['More than one ' modality_data{mod} ' merged file for ' participant{p}]);
                end

                jobfile = {[subject_file_path '/Append_planar_job.m']};
                jobs = repmat(jobfile, 1, 1);
                inputs = cell(0, 1);
                inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]);
                eval(['cd ' outcome_path]) % To ensure it saves them there
                spm('defaults', 'EEG');
                spm_jobman('run', jobs, inputs{:});
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Adjust names
                folders = dir([outcome_path '/']);
                infolder = find(endsWith({folders.name},[modality_data{mod} '_cleaned.mat']) & (...
                            startsWith({folders.name},['P' participant{p}])));
                if length(infolder) > 1
                    error(['More than one cleaned file for ' participant{p}]);
                end      

                % Should only be one, but just in case
                load([outcome_path '/' folders(infolder).name]);
                D.data.fname = [outcome_path '/' participant{p} '_' modality_data{mod} '_cleaned.dat'];
                D.fname = [participant{p} '_' modality_data{mod} '_cleaned.mat'];

                % Adjust name of .dat file with data as well
                datfilename = [outcome_path '/P' participant{p} '_' modality_data{mod} '_cleaned.dat'];  %#ok<*SAGROW> 
                % If it was already renamed there will be not imported_dat
                if exist(datfilename,'file')
                    % Delete original .dat file (as new one will have same name)
                    delete([outcome_path '/' participant{p} '_' modality_data{mod} '_cleaned.dat']);
                    datnewname = [outcome_path '/' participant{p} '_' modality_data{mod} '_cleaned.dat'];  %#ok<*SAGROW>  
                    movefile(datfilename,datnewname);
                end
                
                % Now, save .mat (will overwrite previous .mat)
                save([outcome_path '/' participant{p} '_' modality_data{mod} '_cleaned'],'D');                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                % No need to adjust names (resulting file is named the same)
                disp(' ');      
                disp('-------------------------');
                disp([participant{p} ' COMBINED PLANAR ' modality_data{mod}]);
                disp(datetime)
                disp(' ');                 
            end
        end

        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_cropch';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 
        
        % Do not delete previous steps (cleaned)

    end   
end

%% Crop channels that won't be used in each modality

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

%% Average across conditions (DEV, STD, MMN)

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

%% Obtain Mismatch

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
            D.data.fname = [outcome_path '/' participant{p} '_' modality_data{mod} '_MMN.dat'];
            D.fname = [participant{p} '_' modality_data{mod} '_MMN.mat'];

            % Adjust name of .dat file with data as well
            datfilename = [outcome_path '/W' participant{p} '_' modality_data{mod} '_averaged.dat'];  %#ok<*SAGROW> 
            % If it was already renamed there will be not imported_dat
            if exist(datfilename,'file')
                datnewname = [outcome_path '/' participant{p} '_' modality_data{mod} '_MMN.dat'];  %#ok<*SAGROW>  
                movefile(datfilename,datnewname);
            end

            % Now, save .mat (will overwrite previous .mat)
            save([outcome_path '/' participant{p} '_' modality_data{mod} '_MMN'],'D');   
            % Delete original file
            delete([outcome_path '/W' participant{p} '_' modality_data{mod} '_averaged.mat']);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % No need to adjust names (resulting file is named the same)
            disp(' ');      
            disp('-------------------------');
            disp(['MMN file computed for ' participant{p} ' ' modality_data{mod}]);
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

%% Compute forward model

% To take into account in SPM forward models (as in manual):
% 1) ~10,000 vertices per hemisphere using the highest (fine) resolution
% 2) The individual cortical mesh is obtained automatically from a canonical mesh in MNI
% 3) Individual MRIs are not warped to MNI space, only the mesh with source
% vertices is inverse normalized from an MNI template to the native MRI,
% which is why, upon averaging sources, number of vertices match and
% average is possible (as well as statistics)

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
            % Originally we said 'mni', but after realizing the MRI is not in MNI space yet, we go for 'world'
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
            infolder = find(endsWith({folders.name},[modality_data{mod} '_MMN.mat']) & (...
                        startsWith({folders.name},[participant{p}])));
            if isempty(infolder)
                error(['No ' modality_data{mod} ' MMN files for ' participant{p}]);
            end

            % Should not be by definition, as it is 
            if length(infolder) > 1
                error(['More than one ' modality_data{mod} ' MMN file for ' participant{p}]);
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
        subject_array{pos_subj,3} = 'needs_inverse';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 

    end   
end

%% (Manual) revise corregistrations

%% Compute inverse solutions

% To take into account in SPM inverse solutions (as in manual):
% * Based on an empirical Bayesian formalism, the inversion is meant to be generic in the sense
% it can incorporate and estimate the relevance of multiple constraints of varied nature; 
% datadriven relevance estimation being made possible through Bayesian model comparison
% In other words: the inverse reconstruction step consists in Bayesian inversion

% Individual MRIs are not warped to MNI space, only the mesh with source
% vertices is inverse normalized from an MNI template to the native MRI,
% which is why, upon averaging sources, number of vertices match and
% average is possible (as well as statistics). All of this was tested.

% In SPM we can use two types of inverse solution algorithms as per manual:
% 1) Minimum Norm, but called Independent and Identically-Distributed (IID) in SPM (IID is what they consider the prior)
% 2) Multiple Sparse Priors (MSP), unique to SPM, which corresponds to a sparse prior on the sources, namely that
% only a few are active. Specifically, we use Greedy Search (GS), one of several fitting algorithms for optimizing 
% the MSP approach. GS is the quickest.

% For EEG file: 1) EEG IID 2) EEG GS
% For MEG file: 1) MEG and PLANAR IID 2) MEG IID 3) PLANAR IID 4) MEG and PLANAR GS 5) MEG GS 6) PLANAR GS
% For BIMODAL file: 1) EEG and MEG and PLANAR IID 2) EEG and MEG IID 3) EEG and PLANAR IID 
% 4) EEG and MEG and PLANAR GS 5) EEG and MEG GS 6) EEG and PLANAR GS

for p = 1:length(participant)
    % Check log info about the subject
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(subject_array{pos_subj,3},'needs_inverse')
        
        % For every new participant (MMN file) stablish order of display
         inv_sol_order_eeg = 1;
         inv_sol_order_meg = 1;
         inv_sol_order_bimodal = 1;
        
        % Do for every inverse solution algorithm
        for isa = 1:length(Inverse_solution_algorithm)
            
            % Compute forward model for all modalities
            for mod = 1:length(modality_data)

                
                % There will be no EEG file to obtain MMN from
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

                % Determine combinations of sensors based on modality
                if strcmp(modality_data{mod},'EEG')
                    for esc = 1:length(EEG_sensor_combinations)
                        jobfile = {[subject_file_path '/Inverse_solution_' EEG_sensor_combinations{esc} '_' Inverse_solution_algorithm{isa} '_job.m']};
                        jobs = repmat(jobfile, 1, 1);
                        inputs = cell(0, 1);
                        inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]);
                        inputs{1,2} = inv_sol_order_eeg;
                        eval(['cd ' outcome_path]) % To ensure it saves files there
                        spm('defaults', 'EEG');
                        spm_jobman('run', jobs, inputs{:});
                        inv_sol_order_eeg = inv_sol_order_eeg + 1;

                        % No need to adjust names (resulting file is named the same)
                        disp(' ');      
                        disp('-------------------------');
                        disp(['Inverse modeling computed for ' participant{p} ' ' EEG_sensor_combinations{esc} ' ' Inverse_solution_algorithm{isa}]);
                        disp(datetime)
                        disp(' ');
                    end
                elseif strcmp(modality_data{mod},'MEG')
                    for esc = 1:length(MEG_sensor_combinations)
                        jobfile = {[subject_file_path '/Inverse_solution_' MEG_sensor_combinations{esc} '_' Inverse_solution_algorithm{isa} '_job.m']};
                        jobs = repmat(jobfile, 1, 1);
                        inputs = cell(0, 1);
                        inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]);
                        inputs{1,2} = inv_sol_order_meg;
                        eval(['cd ' outcome_path]) % To ensure it saves files there
                        spm('defaults', 'EEG');
                        spm_jobman('run', jobs, inputs{:});
                        inv_sol_order_meg = inv_sol_order_meg + 1;

                        % No need to adjust names (resulting file is named the same)
                        disp(' ');      
                        disp('-------------------------');
                        disp(['Inverse modeling computed for ' participant{p} ' ' MEG_sensor_combinations{esc} ' ' Inverse_solution_algorithm{isa}]);
                        disp(datetime)
                        disp(' ');
                    end
                elseif strcmp(modality_data{mod},'BIMODAL')
                    for esc = 1:length(BIMODAL_sensor_combinations)
                        jobfile = {[subject_file_path '/Inverse_solution_' BIMODAL_sensor_combinations{esc} '_' Inverse_solution_algorithm{isa} '_job.m']};
                        jobs = repmat(jobfile, 1, 1);
                        inputs = cell(0, 1);
                        inputs{1,1} = cellstr([outcome_path '/' folders(infolder).name]);
                        inputs{1,2} = inv_sol_order_bimodal;
                        eval(['cd ' outcome_path]) % To ensure it saves files there
                        spm('defaults', 'EEG');
                        spm_jobman('run', jobs, inputs{:});
                        inv_sol_order_bimodal = inv_sol_order_bimodal + 1;

                        % No need to adjust names (resulting file is named the same)
                        disp(' ');      
                        disp('-------------------------');
                        disp(['Inverse modeling computed for ' participant{p} ' ' BIMODAL_sensor_combinations{esc} ' ' Inverse_solution_algorithm{isa}]);
                        disp(datetime)
                        disp(' ');
                    end
                end

                % No need to delete any file or change names
                % (as it will include inverse model in same file)
            end
        end
        % If successful, update subject_array for this subject
        subject_array{pos_subj,3} = 'needs_grandaverage';
        save([subject_file_path '/' subject_array_filename '.mat'],'subject_array') 

    end   
end

%% Average of scalp and sources in SPM

if strcmp(compute_grand_average, 'YES')
for pg = 1:length(participant_group)
        
    inputs_eeg = {};
    inputs_meg = {};
    inputs_bimodal = {};
    
    for p = 1:length(participant)

        % Check log info about the subject
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        if ~strcmp(subject_array{pos_subj,3},'needs_grandaverage')
            continue
        end
        
        % Include only participants that correspond to the group
        if ~strcmp(subject_array{pos_subj,2},participant_group{pg})
           continue; 
        end

        % Compute forward model for all modalities
        for mod = 1:length(modality_data)

            % There will be no EEG file to obtain MMN from
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

            % Store in inputs variable
            if strcmp(modality_data{mod},'EEG')
                inputs_eeg{length(inputs_eeg)+1,1} = [outcome_path '/' folders(infolder).name];
            elseif strcmp(modality_data{mod},'MEG')
                inputs_meg{length(inputs_meg)+1,1} = [outcome_path '/' folders(infolder).name];
            elseif strcmp(modality_data{mod},'BIMODAL')
                inputs_bimodal{length(inputs_bimodal)+1,1} = [outcome_path '/' folders(infolder).name];
            end
        end
    end
    
    has_eeg = find(contains(modality_data,'EEG'));
    
    if ~isempty(has_eeg) % If this was one of the choices
        
        disp(' ');      
        disp('-------------------------');
        disp(['COMPUTING EEG GRAND AVERAGE(' participant_group{pg} ')']);
        disp(datetime)
        disp(' '); 
        
        % Test to be absolutely sure about sensor correspondence
        load('/private/path/PEPP/User/DCM/Scripts_final_pipeline/list_EEG_sensors.mat')
        for i = 1:length(inputs_eeg)
            load(inputs_eeg{i});
            if contains(inputs_eeg{i},'2218A_EEG') 
                try % Will fail if pos 64 was already corrected
                    D.sensors.eeg.chanpos(64,:) = [];
                    D.sensors.eeg.chantype(64) = [];
                    D.sensors.eeg.chanunit(64) = [];
                    D.sensors.eeg.elecpos(64,:) = [];
                    D.sensors.eeg.label(64) = [];
                    D.sensors.eeg.tra(64,:) = [];
                    D.sensors.eeg.tra(:,64) = [];
                    D.sensors.eeg.label{61} = 'EEG061';
                    D.sensors.eeg.label{62} = 'EEG063';
                    D.sensors.eeg.label{63} = 'EEG064';
                    save(inputs_eeg{i},'D');
                catch
                    disp(' ');      
                    disp('-------------------------');
                    disp('2218A EEG sensors already corrected');
                    disp(datetime)
                    disp(' '); 
                end
            end
            disp(inputs_eeg{i});
            if length(D.channels) ~=60; error('here');end
            if length(find(ismember({D.channels.label},list_EEG_sensors))) ~=60; error('here');end
            if ~isfield(D.sensors,'eeg'); error('here'); end
            if length(D.sensors.eeg.label) ~= 63; error('here'); end
            if ~isfield(D.sensors,'meg'); error('here');end
            if length(D.sensors.meg.label) ~= 306; error('here');end
        end
        
        % EEG average
        jobfile = {[subject_file_path '/Grandmean_EEG_' participant_group{pg} '_job.m']};
        jobs = repmat(jobfile, 1, 1);
        spm('defaults', 'EEG');
        spm_jobman('run', jobs, inputs_eeg);
        
    end
    
    has_meg = find(contains(modality_data,'MEG'));
    
    if ~isempty(has_meg) % If this was one of the choices
        
        disp(' ');      
        disp('-------------------------');
        disp(['COMPUTING MEG GRAND AVERAGE(' participant_group{pg} ')']);
        disp(datetime)
        disp(' '); 
        
        % For MEG cropped files only, if they have D.sensor.eeg field it
        % won't merge them with subjects that do not have it (subjects
        % without EEG chans), so it is best to remove it from everyone
        % before
        for i = 1:length(inputs_meg)
            load(inputs_meg{i});
            if isfield(D.sensors,'eeg')
                D.sensors = rmfield(D.sensors,'eeg');
                save(inputs_meg{i},'D');
            end
        end
        
        % MEG average
        jobfile = {[subject_file_path '/Grandmean_MEG_' participant_group{pg} '_job.m']};
        jobs = repmat(jobfile, 1, 1);
        spm('defaults', 'EEG');
        spm_jobman('run', jobs, inputs_meg);    
        
        
    end
    
    has_bimodal = find(contains(modality_data,'BIMODAL'));
    
    if ~isempty(has_bimodal) % If this was one of the choices
        
        disp(' ');      
        disp('-------------------------');
        disp(['COMPUTING BIMODAL GRAND AVERAGE(' participant_group{pg} ')']);
        disp(datetime)
        disp(' '); 
        
        % Fix 2218A
        for i = 1:length(inputs_bimodal)
            if contains(inputs_bimodal{i},'2218A_BIMODAL') 
                load(inputs_bimodal{i});
                try % Will fail if pos 64 was already corrected
                    D.sensors.eeg.chanpos(64,:) = [];
                    D.sensors.eeg.chantype(64) = [];
                    D.sensors.eeg.chanunit(64) = [];
                    D.sensors.eeg.elecpos(64,:) = [];
                    D.sensors.eeg.label(64) = [];
                    D.sensors.eeg.tra(64,:) = [];
                    D.sensors.eeg.tra(:,64) = [];
                    D.sensors.eeg.label{61} = 'EEG061';
                    D.sensors.eeg.label{62} = 'EEG063';
                    D.sensors.eeg.label{63} = 'EEG064';
                    save(inputs_bimodal{i},'D');
                catch
                    disp(' ');      
                    disp('-------------------------');
                    disp('2218A bimodal already corrected');
                    disp(datetime)
                    disp(' '); 
                end
            end
        end
        
        
        % BIMODAL average
        jobfile = {[subject_file_path '/Grandmean_BIMODAL_' participant_group{pg} '_job.m']};
        jobs = repmat(jobfile, 1, 1);
        spm('defaults', 'EEG');
        spm_jobman('run', jobs, inputs_bimodal);

    end
    
    % No need to delete any file or change names
    % (as it will create a new one)
end
end

%% Ready for DCM
