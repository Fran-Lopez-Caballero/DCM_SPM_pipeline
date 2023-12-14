%% Compute Bayesian Model Reduction (NMDA)

%% Define variables
clear; clc;
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
load([subject_file_path '/' subject_array_filename '.mat'],'subject_array');
modality_data = 'EEG'; % Choose one at a time ('EEG' or 'MEG'
MMN_types = 'pMMN'; % Choose one at a time (pMMN or dMMN)
MMN_combinations_DCM = {[1,2],[1,3]}; % 1 = STD, 2 = PDev, 3 = DDev; so [1,2] if pMMN and [1,3] if dMMN
if strcmp(MMN_types, 'pMMN')
    mmnty = 1; % For next line
elseif strcmp(MMN_types,'dMMN')
    mmnty = 2; % For next line
end
DCM_type = 'NMDA'; % 'CMC' or 'NMDA'
tag_DCM_solved = 'solved_v3'; % 'solved_v3' or '_link'
CRITERIA_R2 = 85; % participants with fits under this R2 will be attempted to be re-inversed
% Define ROIs
rois = {'left_A1','right_A1','left_STG','right_STG','left_IFG','right_IFG'};
mni_left_A1 = [-42, -22,7];
mni_right_A1 =  [46, -14, 8];
mni_left_STG = [-61, -32, 8];
mni_right_STG = [59, -25, 8];
mni_left_IFG = [-46, 20, 8];
mni_right_IFG = [46, 20, 8];

% R2=100*Predicted response/(Predicted response + Residual response). Response being PCA first eight components
participant = {subject_array{:,1}}; %#ok<*CCAT1,*USENS>
iteration = 2; % Number of times you have run this script
spm('defaults', 'EEG');

%% Load structure with non-inverted DCM for each subject (most complex model defined by default)

% This step is already done in A2_Prepare_DCM_struc, so just load the file

load([outcome_path '/DCM_structs/DCMs_' modality_data '_' MMN_types '_' DCM_type '.mat'],'DCMs'); 
load([outcome_path '/DCM_structs/Participant_index_' modality_data '_' MMN_types '_' DCM_type '.mat'],'participant_index');
DCMs_not_solved = DCMs;

%% Obtain structure with inverted DCM for each subject (most complex model)

DCMs = {};
participant_index_2 = {};
participant_group = {};
% Load solved 
for p = 1:length(participant)
    
    % Load DCM file for that subject
    try
        load([outcome_path '/' participant{p} '_' modality_data '_DCM_' DCM_type '_' MMN_types '_' tag_DCM_solved '.mat'],'DCM');
    catch
        error(['No ' participant{p} '_' modality_data '_DCM_' DCM_type '_' MMN_types '_' tag_DCM_solved '.mat file found']);
    end

    DCMs{p,1} = DCM; %#ok<*NODEF,*AGROW>
    participant_index_2{p,1} = participant{p}; %#ok<*SAGROW>
    participant_index_2{p,1} = participant{p}; %#ok<*SAGROW>
    participant_group{p,1} = subject_array{p,2};
end

% Eliminate empty cells (common issue with EEG/MEG alternative)
DCMs = DCMs(~cellfun('isempty',DCMs));
participant_index_2 = participant_index_2(~cellfun('isempty',participant_index_2));

if ~isequal(participant_index,participant_index_2)
    error('subject arrays are not equal between solved and non-solved DCMs')
end

if length(participant_index_2) ~= length(DCMs)
    error('participant list does not match DCM list');
end

%% free_parameters_besides_B = {'T'}; 

% For here on out, script is based on template from Ryszard Auksztulewicz

free_parameters_besides_B = {'T' 'BN{1}' 'A{1}' 'A{2}' 'AN{1}' 'AN{2}'}; 
% A matrix describing connections that are constant (common to all conditions, std and dev)
% B differences between conditions (deviant and standard) % MMN
% 'T' 'BN{1}' 'A{1}' 'A{2}' 'AN{1}' 'AN{2}'
% free_parameters_besides_B = {'T'}; 
% allowing more parameters to be optimized at the group level we are
% changing some fits at the T parameter
% adding more degrees of freedom
p_thresh = .99; % probability threshold for plotting; 
% you can make it higher (e.g. 0.99 or 0.995)

%% Debugging BMR

model_space = fliplr(gen_model_space_NMDA());
GCM = cell(length(DCMs),length(model_space));
for k = 1:length(DCMs)
    for i=1:length(model_space)
        if i==1
            GCM{k,i}=DCMs{k}; % Full model
            GCM{k,i}.M.pC = spm_unvec(GCM{k,i}.M.pC,GCM{k,i}.M.pE);
            GCM{k,i}.M.pC.B = GCM{k,i}.B;
            ns = size(GCM{k,i}.M.pC.B{1},1);
            GCM{k,i}.M.pC.B{1}(find(eye(ns)==1))=0;
            for j=1:length(free_parameters_besides_B)
                parname = free_parameters_besides_B{j};
                eval(strcat(['GCM{k,i}.M.pC.',parname,' = GCM{k,i}.M.pC.',parname,'+1;']))
                eval(strcat(['GCM{k,i}.M.pC.',parname,'(find(GCM{k,i}.M.pE.',parname,'==-32))=0;']))
            end
        else
            GCM{k,i}=DCMs_not_solved{k};
            GCM{k,i}.B{1} = model_space{i}.matrix; 
            GCM{k,i}.M = GCM{k,1}.M; 
            GCM{k,i}.M.pC = spm_unvec(GCM{k,i}.M.pC,GCM{k,i}.M.pE);
            GCM{k,i}.M.pC.B{1} = model_space{i}.matrix; % Ryszard
            for j=1:length(free_parameters_besides_B)
                parname = free_parameters_besides_B{j};
                eval(strcat(['GCM{k,i}.M.pC.',parname,' = GCM{k,i}.M.pC.',parname,'+1;']))
                eval(strcat(['GCM{k,i}.M.pC.',parname,'(find(GCM{k,i}.M.pE.',parname,'==-32))=0;']))
            end
        end
        names{i}=model_space{i}.name;
    end
end
[RCM,BMC,BMA] = spm_dcm_bmr(GCM); %#ok<*ASGLU>

% Plot model probabilities by subject
for i=1:length(DCMs)
    model_prob(i,:)=BMC(i).P'; %#ok<*SAGROW>
end
figure();
colormap(flipud(gray));
imagesc(model_prob)
colorbar;
xticklabels(names);
title('Model Probability per Subject')
xlabel('Model')
ylabel('Subject')

%% new code: PEB

no_subjects = size(model_prob,1);
M = [];
M.X = [ones(no_subjects,1) zeros(no_subjects,1)];
% M.X = [ones(no_subjects,1) repelem(-1,no_subjects)'];
M.X(27:end,2) = 1; % patients
% M.X(find(strcmp(participant_group,'FE')),2) = 1; % patients


% the first column contains commonalities in parameters (i.e., parameter
% average across patients and controls); the second colum contains
% differences between patients and controls

clear PEB F_PEB

for i=1:length(model_space)
    [PEB{i},DCM] = spm_dcm_peb(RCM(:,i),M,['B' free_parameters_besides_B]); % PEB contains group models
end

%% run Bayesian model averaging

BMA = spm_dcm_bma(PEB); % average across all models

% convert parameter variance to probability
BMA.Pp = BMA.Cp*0; 

M = BMA.Ep;
S = BMA.Cp;

N = 10000;

for i=1:size(M,1)
    for j=1:size(M,2)
        m = M(i,j);
        s = S(i,j);
        x = s*randn(N,1)+m;
        b = exp(x);
        pdecrease = length(find(b<1))/N;
        pincrease = length(find(b>1))/N;
        BMA.Pp(i,j)=max([pdecrease pincrease]);
    end
end

%% plot

model_names = cell(0);
for i=1:length(model_space)
    model_names{i} = model_space{i}.name;
end

figure;
subplot(2,2,1)
rel_lme = BMA.F - min(BMA.F);
bar(rel_lme,'facecolor',[.5 .5 .5],'edgecolor','none')
hold on
bar(BMA.Mocc,rel_lme(BMA.Mocc),'facecolor',[1 .5 .5],'edgecolor','none')
xlabel('model')
xticks(1:length(model_space))
xticklabels(model_names)
ylabel('relative log-model evidence')
title('model comparison')
legend({'' 'models for BMA'})

winning_model = find(BMA.F==max(BMA.F));

subplot(2,2,2)
hold on
bar(BMA.Ep(:,1),'facecolor',[.5 .5 .5],'edgecolor','none')
errorbar(BMA.Ep(:,1),BMA.Cp(:,1),'linestyle','none','linewidth',2,'capsize',0)
axisrange = max(abs(BMA.Ep(:,1))+abs(BMA.Cp(:,1)))*1.5;
scatter(find(BMA.Pp(:,1)>p_thresh),repelem(axisrange*.9,length(find(BMA.Pp(:,1)>p_thresh))),15,'k*')
% scatter(find(BMA.Pp(:,1)>p_thresh),axisrange*.9,15,'k*')
xticks(1:length(BMA.Ep(:,1)))
xticklabels(PEB{winning_model}.Pnames)
xtickangle(45)
ylim([-axisrange axisrange])
title('group commonalities')

subplot(2,2,4)
hold on
bar(BMA.Ep(:,2),'facecolor',[.5 .5 .5],'edgecolor','none')
errorbar(BMA.Ep(:,2),BMA.Cp(:,2),'linestyle','none','linewidth',2,'capsize',0)
axisrange = max(abs(BMA.Ep(:,2))+abs(BMA.Cp(:,2)))*1.5;
scatter(find(BMA.Pp(:,2)>p_thresh),repelem(axisrange*.9,length(find(BMA.Pp(:,2)>p_thresh))),15,'k*')
% scatter(find(BMA.Pp(:,2)>p_thresh),axisrange*.9,15,'k*')
xticks(1:length(BMA.Ep(:,2)))
xticklabels(PEB{winning_model}.Pnames)
xtickangle(45)
ylim([-axisrange axisrange])
title('group differences')

%% Results with gui

% spm_dcm_peb_review(PEB{1},RCM)

