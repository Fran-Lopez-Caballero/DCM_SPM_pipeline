%% Compute R2 (square of correlation) between predicted and real data for each DCM

%% Define variables and start SPM

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

set(0,'defaultfigurecolor',[1 1 1]); % I want white backgrounds
% participant_group = {'C','FE'};
subject_array_filename = 'subject_array_EEG_matched';
load([subject_file_path '/' subject_array_filename '.mat']);
modality_data = {'EEG'}; % 'EEG','MEG'
MMN_types = {'pMMN'}; % 'pMMN','dMMN'
participant = {subject_array{:,1}};
save_figures = 'NO'; % 'YES' or 'NO'
tag_DCM_solved = '_v3';

%% Compute R2 and plot an histogram

for mode = 1:length(modality_data)
for tymmn = 1:length(MMN_types)
    R2_single = [];
    good_sub_pos = [];gsp_pos = 1;
    
    for p = 1:length(participant)

        % Choose based on whether there will be files or not
        pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
        if strcmp(modality_data{mode},'EEG') && (strcmp(subject_array{pos_subj,4},'no_EEG_chans') || strcmp(subject_array{pos_subj,5},'bad_EEG'))
            continue;
        end
        if strcmp(modality_data{mode},'MEG') && strcmp(subject_array{pos_subj,6},'bad_MEG')
            continue;
        end

        % Adjustement for 2379 not having pMMN MEG DCM
        if strcmp(modality_data{mode},'MEG') && strcmp(MMN_types{tymmn},'pMMN') && strcmp(subject_array{pos_subj,14},'no_pMMN_MEG_CMC_DCM')
            continue;
        end
        
        % Load file from that subject
        load([outcome_path '/' participant{p} '_' modality_data{mode} '_DCM_CMC_' MMN_types{tymmn} '_solved' tag_DCM_solved '.mat']);
        
        RSS=sum(sum([DCM.R{1}.^2 DCM.R{2}.^2])); % Residual?
        PSS=sum(sum([(DCM.H{1}+DCM.R{1}).^2 (DCM.H{2}+DCM.R{2}).^2])); % Predicted?
        R2_single(p)=100*PSS/(PSS + RSS);
        if R2_single(p) < 85
            % disp([participant{p} ' '  modality_data{mode} ' ' MMN_types{tymmn} ' has ' num2str(R2_single(p)) ' fit'])
            disp(participant{p});
            if strcmp(MMN_types{tymmn},'pMMN')
                subject_array{p,21} = 'needs_CMC_rerun_pMMN';
            elseif strcmp(MMN_types{tymmn},'dMMN')
                subject_array{p,22} = 'needs_CMC_rerun_dMMN';
            end
        end
        good_sub_pos(gsp_pos) = p;
        gsp_pos = gsp_pos +1;
        
        if any(isnan(spm_vec(DCM.Ep))) || any(isinf(spm_vec(DCM.Ep))) || any(isnan(spm_vec(DCM.Cp))) || any(isinf(spm_vec(DCM.Cp)))
            nans(p)=1;
        end
    end
    
    figure()
    R2_single = R2_single(good_sub_pos);
    participant_header = participant(good_sub_pos);
    display_R2_single = {};
    if length(participant_header) ~= length(R2_single)
        error('dimension mismatch between header and R2 results');
    end
    for rwo = 1:length(participant_header)
        display_R2_single{rwo,1} = participant_header(rwo);
        display_R2_single{rwo,2} = R2_single(rwo); %#ok<*SAGROW>
    end
    eval([MMN_types{tymmn} '_R2_single = display_R2_single;'])
    histogram(R2_single,30,'FaceColor','blue');
    ylim([0 10])
    title(['CMC ' modality_data{mode} ' ' MMN_types{tymmn} ' R^2 Individual DCMs']);
    if strcmp(save_figures,'YES')
        saveas(gcf,[outcome_path 'Fits_DCM/histogram_CMC_' modality_data{mode} '_' MMN_types{tymmn} '.png']);
        close all;
    end
    R2_single_avg = mean(R2_single);
    disp(['Average ' MMN_types{tymmn} ' ' modality_data{mode} ' fit is ' num2str(R2_single_avg)]);
end
end

%% Plot first principal Eigenvalue (fit between predicted and real data) and save figures

cols_unsort  = {'red','blue','black','green'};       
count=0;

for mode = 1:length(modality_data)
for tymmn = 1:length(MMN_types)
for p = 1:length(participant)
    % Choose based on whether there will be files or not
    pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
    if strcmp(modality_data{mode},'EEG') && (strcmp(subject_array{pos_subj,4},'no_EEG_chans') || strcmp(subject_array{pos_subj,5},'bad_EEG'))
        continue;
    end
    if strcmp(modality_data{mode},'MEG') && strcmp(subject_array{pos_subj,6},'bad_MEG')
        continue;
    end

    if strcmp(modality_data{mode},'MEG') && strcmp(MMN_types{tymmn},'pMMN') && strcmp(subject_array{pos_subj,14},'no_pMMN_MEG_CMC_DCM')
        continue;
    end
    
    % Load DCM file for that subject
    try
        load([outcome_path '/' participant{p} '_' modality_data{mode} '_DCM_CMC_' MMN_types{tymmn} '_solved' tag_DCM_solved '.mat']);
    catch
        error(['No ' participant{p} '_' modality_data{mode} '_DCM_CMC_' MMN_types{tymmn} '_solved' tag_DCM_solved '.mat file found']);
    end
    
    figure;
    for c = 1:2
        plot(DCM.H{c}(:,1), 'color', cols_unsort{c}, 'Linewidth', 1.5, 'LineStyle','--'); hold on
        plot(DCM.H{c}(:,1) + DCM.R{c}(:,1), 'color', cols_unsort{c}, 'Linewidth', 1.5);
        ylim([-10 10]);           
        set(gcf, 'color', 'w');
    end
    title(strcat([modality_data{mode} ' CMC ' MMN_types{tymmn} ' ' participant{p}]),'Interpreter', 'none')
    legend({'Std(pred)', 'Std(obs)', 'Dev(pred)', 'Dev(obs)'});
    Image = getframe(gcf);    
    if strcmp(save_figures,'YES')
        saveas(gcf,[outcome_path 'Fits_DCM/' modality_data{mode} '_DCM_CMC_' MMN_types{tymmn} '_' participant{p} '.png']);
        close all;
    end
end
end
end
