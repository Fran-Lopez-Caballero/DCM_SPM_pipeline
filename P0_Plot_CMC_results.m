%% Translate DCM structure results into meaningful variables (to interpret)
% General advice: when I write MMN here it refers to the conditional change in 
% connections from baseline(standard) to deviant. If willing to obtain the 
% connectivity strength to the deviant alone, one would have to sum A + B or M + N

%% Define variables and start SPM

clear;clc;

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
subject_array_filename = 'subject_array_EEG_matched';
load([subject_file_path '/' subject_array_filename '.mat']);
modality_data = {'EEG'}; % {'EEG','MEG'};
MMN_types = {'dMMN'}; % {'pMMN','dMMN'};
ROIs = {'A1_L', 'A1_R', 'STG_L', 'STG_R', 'IFG_L', 'IFG_R'};
participant = {subject_array{:,1}};
intrinsic_connections = {'sp2sp', 'sp2ss', 'ii2ss'};
extract_variables = 'YES'; % 'YES' or 'NO'
plot_butterfly = 'NO';
plot_gavr = 'NO';
plot_connectivity = 'YES';
tag_DCM_solved = '_v3';

%% Loop through subjects, modality of data and mmn type

if strcmp(extract_variables,'YES')
for mode = 1:length(modality_data)
for tymmn = 1:length(MMN_types)
    
PEPP_MMN_Baseline_DCM = {};
% Create storing variable
if 1 == 1 % Silly way to compress section
PEPP_MMN_Baseline_DCM{1,1} = 'MEASURE/PARTICIPANT';
PEPP_MMN_Baseline_DCM{3,1} = '----SOURCES----';
PEPP_MMN_Baseline_DCM{4,1} = 'STD_A1_L_ii';
PEPP_MMN_Baseline_DCM{5,1} = 'STD_A1_R_ii';
PEPP_MMN_Baseline_DCM{6,1} = 'STD_STG_L_ii';
PEPP_MMN_Baseline_DCM{7,1} = 'STD_STG_R_ii';
PEPP_MMN_Baseline_DCM{8,1} = 'STD_IFG_L_ii';
PEPP_MMN_Baseline_DCM{9,1} = 'STD_IFG_R_ii';
PEPP_MMN_Baseline_DCM{10,1} = 'STD_A1_L_ei';
PEPP_MMN_Baseline_DCM{11,1} = 'STD_A1_R_ei';
PEPP_MMN_Baseline_DCM{12,1} = 'STD_STG_L_ei';
PEPP_MMN_Baseline_DCM{13,1} = 'STD_STG_R_ei';
PEPP_MMN_Baseline_DCM{14,1} = 'STD_IFG_L_ei';
PEPP_MMN_Baseline_DCM{15,1} = 'STD_IFG_R_ei';
PEPP_MMN_Baseline_DCM{16,1} = 'STD_A1_L_py';
PEPP_MMN_Baseline_DCM{17,1} = 'STD_A1_R_py';
PEPP_MMN_Baseline_DCM{18,1} = 'STD_STG_L_py';
PEPP_MMN_Baseline_DCM{19,1} = 'STD_STG_R_py';
PEPP_MMN_Baseline_DCM{20,1} = 'STD_IFG_L_py';
PEPP_MMN_Baseline_DCM{21,1} = 'STD_IFG_R_py';
PEPP_MMN_Baseline_DCM{22,1} = 'DEV_A1_L_ii';
PEPP_MMN_Baseline_DCM{23,1} = 'DEV_A1_R_ii';
PEPP_MMN_Baseline_DCM{24,1} = 'DEV_STG_L_ii';
PEPP_MMN_Baseline_DCM{25,1} = 'DEV_STG_R_ii';
PEPP_MMN_Baseline_DCM{26,1} = 'DEV_IFG_L_ii';
PEPP_MMN_Baseline_DCM{27,1} = 'DEV_IFG_R_ii';
PEPP_MMN_Baseline_DCM{28,1} = 'DEV_A1_L_ei';
PEPP_MMN_Baseline_DCM{29,1} = 'DEV_A1_R_ei';
PEPP_MMN_Baseline_DCM{30,1} = 'DEV_STG_L_ei';
PEPP_MMN_Baseline_DCM{31,1} = 'DEV_STG_R_ei';
PEPP_MMN_Baseline_DCM{32,1} = 'DEV_IFG_L_ei';
PEPP_MMN_Baseline_DCM{33,1} = 'DEV_IFG_R_ei';
PEPP_MMN_Baseline_DCM{34,1} = 'DEV_A1_L_py';
PEPP_MMN_Baseline_DCM{35,1} = 'DEV_A1_R_py';
PEPP_MMN_Baseline_DCM{36,1} = 'DEV_STG_L_py';
PEPP_MMN_Baseline_DCM{37,1} = 'DEV_STG_R_py';
PEPP_MMN_Baseline_DCM{38,1} = 'DEV_IFG_L_py';
PEPP_MMN_Baseline_DCM{39,1} = 'DEV_IFG_R_py';
PEPP_MMN_Baseline_DCM{40,1} = '----FORWARD CONNECTIVITY----';
PEPP_MMN_Baseline_DCM{41,1} = 'STD_A12STG_L_sp2ss';
PEPP_MMN_Baseline_DCM{42,1} = 'STD_A12STG_R_sp2ss';
PEPP_MMN_Baseline_DCM{43,1} = 'STD_STG2IFG_L_sp2ss';
PEPP_MMN_Baseline_DCM{44,1} = 'STD_STG2IFG_R_sp2ss';
PEPP_MMN_Baseline_DCM{45,1} = 'STD_A12STG_L_sp2dp';
PEPP_MMN_Baseline_DCM{46,1} = 'STD_A12STG_R_sp2dp';
PEPP_MMN_Baseline_DCM{47,1} = 'STD_STG2IFG_L_sp2dp';
PEPP_MMN_Baseline_DCM{48,1} = 'STD_STG2IFG_R_sp2dp';
PEPP_MMN_Baseline_DCM{49,1} = 'DEV_A12STG_L_sp2ss';
PEPP_MMN_Baseline_DCM{50,1} = 'DEV_A12STG_R_sp2ss';
PEPP_MMN_Baseline_DCM{51,1} = 'DEV_STG2IFG_L_sp2ss';
PEPP_MMN_Baseline_DCM{52,1} = 'DEV_STG2IFG_R_sp2ss';
PEPP_MMN_Baseline_DCM{53,1} = 'DEV_A12STG_L_sp2dp';
PEPP_MMN_Baseline_DCM{54,1} = 'DEV_A12STG_R_sp2dp';
PEPP_MMN_Baseline_DCM{55,1} = 'DEV_STG2IFG_L_sp2dp';
PEPP_MMN_Baseline_DCM{56,1} = 'DEV_STG2IFG_R_sp2dp';
PEPP_MMN_Baseline_DCM{57,1} = 'MMN_A12STG_L';
PEPP_MMN_Baseline_DCM{58,1} = 'MMN_A12STG_R';
PEPP_MMN_Baseline_DCM{59,1} = 'MMN_STG2IFG_L';
PEPP_MMN_Baseline_DCM{60,1} = 'MMN_STG2IFG_R';
PEPP_MMN_Baseline_DCM{61,1} = '----BACKWARDS CONNECTIVITY----';
PEPP_MMN_Baseline_DCM{62,1} = 'STD_STG2A1_L_dp2sp';
PEPP_MMN_Baseline_DCM{63,1} = 'STD_STG2A1_R_dp2sp';
PEPP_MMN_Baseline_DCM{64,1} = 'STD_IFG2STG_L_dp2sp';
PEPP_MMN_Baseline_DCM{65,1} = 'STD_IFG2STG_R_dp2sp';
PEPP_MMN_Baseline_DCM{66,1} = 'STD_STG2A1_L_dp2ii';
PEPP_MMN_Baseline_DCM{67,1} = 'STD_STG2A1_R_dp2ii';
PEPP_MMN_Baseline_DCM{68,1} = 'STD_IFG2STG_L_dp2ii';
PEPP_MMN_Baseline_DCM{69,1} = 'STD_IFG2STG_R_dp2ii';
PEPP_MMN_Baseline_DCM{70,1} = 'DEV_STG2A1_L_dp2sp';
PEPP_MMN_Baseline_DCM{71,1} = 'DEV_STG2A1_R_dp2sp';
PEPP_MMN_Baseline_DCM{72,1} = 'DEV_IFG2STG_L_dp2sp';
PEPP_MMN_Baseline_DCM{73,1} = 'DEV_IFG2STG_R_dp2sp';
PEPP_MMN_Baseline_DCM{74,1} = 'DEV_STG2A1_L_dp2ii';
PEPP_MMN_Baseline_DCM{75,1} = 'DEV_STG2A1_R_dp2ii';
PEPP_MMN_Baseline_DCM{76,1} = 'DEV_IFG2STG_L_dp2ii';
PEPP_MMN_Baseline_DCM{77,1} = 'DEV_IFG2STG_R_dp2ii';
PEPP_MMN_Baseline_DCM{78,1} = 'MMN_STG2A1_L';
PEPP_MMN_Baseline_DCM{79,1} = 'MMN_STG2A1_R';
PEPP_MMN_Baseline_DCM{80,1} = 'MMN_IFG2STG_L';
PEPP_MMN_Baseline_DCM{81,1} = 'MMN_IFG2STG_R';
PEPP_MMN_Baseline_DCM{82,1} = '----MODULATORY CONNECTIVITY----';
PEPP_MMN_Baseline_DCM{83,1} = 'STD_A1_L_mod';
PEPP_MMN_Baseline_DCM{84,1} = 'STD_A1_R_mod';
PEPP_MMN_Baseline_DCM{85,1} = 'STD_STG_L_mod';
PEPP_MMN_Baseline_DCM{86,1} = 'STD_STG_R_mod';
PEPP_MMN_Baseline_DCM{87,1} = 'STD_IFG_L_mod';
PEPP_MMN_Baseline_DCM{88,1} = 'STD_IFG_R_mod';
PEPP_MMN_Baseline_DCM{89,1} = 'DEV_A1_L_mod';
PEPP_MMN_Baseline_DCM{90,1} = 'DEV_A1_R_mod';
PEPP_MMN_Baseline_DCM{91,1} = 'DEV_STG_L_mod';
PEPP_MMN_Baseline_DCM{92,1} = 'DEV_STG_R_mod';
PEPP_MMN_Baseline_DCM{93,1} = 'DEV_IFG_L_mod';
PEPP_MMN_Baseline_DCM{94,1} = 'DEV_IFG_R_mod';
PEPP_MMN_Baseline_DCM{95,1} = 'MMN_A1_L_mod';
PEPP_MMN_Baseline_DCM{96,1} = 'MMN_A1_R_mod';
PEPP_MMN_Baseline_DCM{97,1} = 'MMN_STG_L_mod';
PEPP_MMN_Baseline_DCM{98,1} = 'MMN_STG_R_mod';
PEPP_MMN_Baseline_DCM{99,1} = 'MMN_IFG_L_mod';
PEPP_MMN_Baseline_DCM{100,1} = 'MMN_IFG_R_mod';
PEPP_MMN_Baseline_DCM{101,1} = '----INSTRINSIC CONNECTIONS----';
PEPP_MMN_Baseline_DCM{102,1} = 'CMC_A1_L_ii2ss';
PEPP_MMN_Baseline_DCM{103,1} = 'CMC_A1_L_sp2sp';
PEPP_MMN_Baseline_DCM{104,1} = 'CMC_A1_L_sp2ss';
PEPP_MMN_Baseline_DCM{105,1} = 'CMC_A1_R_ii2ss';
PEPP_MMN_Baseline_DCM{106,1} = 'CMC_A1_R_sp2sp';
PEPP_MMN_Baseline_DCM{107,1} = 'CMC_A1_R_sp2ss';
PEPP_MMN_Baseline_DCM{108,1} = 'CMC_STG_L_ii2ss';
PEPP_MMN_Baseline_DCM{109,1} = 'CMC_STG_L_sp2sp';
PEPP_MMN_Baseline_DCM{110,1} = 'CMC_STG_L_sp2ss';
PEPP_MMN_Baseline_DCM{111,1} = 'CMC_STG_R_ii2ss';
PEPP_MMN_Baseline_DCM{112,1} = 'CMC_STG_R_sp2sp';
PEPP_MMN_Baseline_DCM{113,1} = 'CMC_STG_R_sp2ss';
PEPP_MMN_Baseline_DCM{114,1} = 'CMC_IFG_L_ii2ss';
PEPP_MMN_Baseline_DCM{115,1} = 'CMC_IFG_L_sp2sp';
PEPP_MMN_Baseline_DCM{116,1} = 'CMC_IFG_L_sp2ss';
PEPP_MMN_Baseline_DCM{117,1} = 'CMC_IFG_R_ii2ss';
PEPP_MMN_Baseline_DCM{118,1} = 'CMC_IFG_R_sp2sp';
PEPP_MMN_Baseline_DCM{119,1} = 'CMC_IFG_R_sp2ss';
end
    
label_c = {}; label_fe = {};
count_c = 1; count_fe = 1;
for p = 1:length(participant)

% Choose based on whether there will be files or not
pos_subj = find(strcmp({subject_array{:,1}},participant{p}));
if strcmp(modality_data{mode},'EEG') && (strcmp(subject_array{pos_subj,4},'no_EEG_chans') || strcmp(subject_array{pos_subj,5},'bad_EEG'))
    continue;
end
if strcmp(modality_data{mode},'MEG') && strcmp(subject_array{pos_subj,6},'bad_MEG')
    continue;
end
if ~strcmp(subject_array{pos_subj,3},'needs_DCM_NMDA') && ~strcmp(subject_array{pos_subj,3},'DCM_computed')
    continue;
end
% Adjustement for 2379 not having pMMN MEG DCM
if strcmp(modality_data{mode},'MEG') && strcmp(MMN_types{tymmn},'pMMN') && strcmp(subject_array{pos_subj,14},'no_pMMN_MEG_CMC_DCM')
    continue;
end

% Load DCM file for that subject
try
    load([outcome_path '/' participant{p} '_' modality_data{mode} '_DCM_CMC_' MMN_types{tymmn} '_solved' tag_DCM_solved '.mat']);
catch
    error(['No ' participant{p} '_' modality_data{mode} '_DCM_CMC_' MMN_types{tymmn} '_solved' tag_DCM_solved '.mat file found']);
end

disp(['Extracting DCM values for ' participant{p} ' (' modality_data{mode} ' ' MMN_types{tymmn} ')'])    

% Associate a position in matrix
PEPP_MMN_Baseline_DCM{1,p+1} = participant{p};
% Save column positions if HC or FE for later
if contains(subject_array(strcmp({subject_array{:,1}},participant{p}),2),'FE')
    pos_fe(count_fe) = p+1; count_fe = count_fe +1;
elseif contains(subject_array(strcmp({subject_array{:,1}},participant{p}),2),'C')
    pos_c(count_c) = p+1; count_c = count_c +1;  
end
    
%% Source waveforms for different populations (Double check which ones are inhibitory and excitatory)

% STANDARDS

% Inhibitory interneurons
STD_A1_L_ii = DCM.K{1,1}(:,1);
STD_A1_R_ii = DCM.K{1,1}(:,2);
STD_STG_L_ii = DCM.K{1,1}(:,3);
STD_STG_R_ii = DCM.K{1,1}(:,4);
STD_IFG_L_ii = DCM.K{1,1}(:,5);
STD_IFG_R_ii = DCM.K{1,1}(:,6);
% Excitatory interneurons
STD_A1_L_ei = DCM.K{1,1}(:,7);
STD_A1_R_ei = DCM.K{1,1}(:,8);
STD_STG_L_ei = DCM.K{1,1}(:,9);
STD_STG_R_ei = DCM.K{1,1}(:,10);
STD_IFG_L_ei = DCM.K{1,1}(:,11);
STD_IFG_R_ei = DCM.K{1,1}(:,12);
% Pyramidal cells (sources)
STD_A1_L_py = DCM.K{1,1}(:,13);
STD_A1_R_py = DCM.K{1,1}(:,14);
STD_STG_L_py = DCM.K{1,1}(:,15);
STD_STG_R_py = DCM.K{1,1}(:,16);
STD_IFG_L_py = DCM.K{1,1}(:,17);
STD_IFG_R_py = DCM.K{1,1}(:,18);

% DEVIANTS
% Inhibitory interneurons
DEV_A1_L_ii = DCM.K{1,2}(:,1);
DEV_A1_R_ii = DCM.K{1,2}(:,2);
DEV_STG_L_ii = DCM.K{1,2}(:,3);
DEV_STG_R_ii = DCM.K{1,2}(:,4);
DEV_IFG_L_ii = DCM.K{1,2}(:,5);
DEV_IFG_R_ii = DCM.K{1,2}(:,6);
% Excitatory interneurons
DEV_A1_L_ei = DCM.K{1,2}(:,7);
DEV_A1_R_ei = DCM.K{1,2}(:,8);
DEV_STG_L_ei = DCM.K{1,2}(:,9);
DEV_STG_R_ei = DCM.K{1,2}(:,10);
DEV_IFG_L_ei = DCM.K{1,2}(:,11);
DEV_IFG_R_ei = DCM.K{1,2}(:,12);
% Pyramidal cells (sources)
DEV_A1_L_py = DCM.K{1,2}(:,13);
DEV_A1_R_py = DCM.K{1,2}(:,14);
DEV_STG_L_py = DCM.K{1,2}(:,15);
DEV_STG_R_py = DCM.K{1,2}(:,16);
DEV_IFG_L_py = DCM.K{1,2}(:,17);
DEV_IFG_R_py = DCM.K{1,2}(:,18);

%% Connectivity measures between regions
% Forward connections baseline (standard) and conditional change (MMN)

% Superior pyramidal to spiny stellate
STD_A12STG_L_sp2ss = DCM.Ep.A{1,1}(3,1);
STD_A12STG_R_sp2ss = DCM.Ep.A{1,1}(4,2);
STD_STG2IFG_L_sp2ss = DCM.Ep.A{1,1}(5,3);
STD_STG2IFG_R_sp2ss = DCM.Ep.A{1,1}(6,4);
% Superior pyramidal to deep pyramidal
STD_A12STG_L_sp2dp = DCM.Ep.A{1,2}(3,1);
STD_A12STG_R_sp2dp = DCM.Ep.A{1,2}(4,2);
STD_STG2IFG_L_sp2dp = DCM.Ep.A{1,2}(5,3);
STD_STG2IFG_R_sp2dp = DCM.Ep.A{1,2}(6,4);

% MMN ones (conditional change): do not have a particular population
MMN_A12STG_L = DCM.Ep.B{1,1}(3,1);
MMN_A12STG_R = DCM.Ep.B{1,1}(4,2);
MMN_STG2IFG_L = DCM.Ep.B{1,1}(5,3);
MMN_STG2IFG_R = DCM.Ep.B{1,1}(6,4);

% DEV (Total response)
% Superior pyramidal to spiny stellate
DEV_A12STG_L_sp2ss = DCM.Ep.A{1,1}(3,1) + MMN_A12STG_L;
DEV_A12STG_R_sp2ss = DCM.Ep.A{1,1}(4,2) + MMN_A12STG_R;
DEV_STG2IFG_L_sp2ss = DCM.Ep.A{1,1}(5,3) + MMN_STG2IFG_L;
DEV_STG2IFG_R_sp2ss = DCM.Ep.A{1,1}(6,4) + MMN_STG2IFG_R;
% Superior pyramidal to deep pyramidal
DEV_A12STG_L_sp2dp = DCM.Ep.A{1,2}(3,1) + MMN_A12STG_L;
DEV_A12STG_R_sp2dp = DCM.Ep.A{1,2}(4,2) + MMN_A12STG_R;
DEV_STG2IFG_L_sp2dp = DCM.Ep.A{1,2}(5,3) + MMN_STG2IFG_L;
DEV_STG2IFG_R_sp2dp = DCM.Ep.A{1,2}(6,4) + MMN_STG2IFG_R;

% Backward connections baseline (standard) and conditional change (MMN)

% Deep pyramidal to superficial pyramidal
STD_STG2A1_L_dp2sp = DCM.Ep.A{1,3}(1,3);
STD_STG2A1_R_dp2sp = DCM.Ep.A{1,3}(2,4);
STD_IFG2STG_L_dp2sp = DCM.Ep.A{1,3}(3,5);
STD_IFG2STG_R_dp2sp = DCM.Ep.A{1,3}(4,6);
% Deep pyramidal to inhibitory interneurons
STD_STG2A1_L_dp2ii = DCM.Ep.A{1,4}(1,3);
STD_STG2A1_R_dp2ii = DCM.Ep.A{1,4}(2,4);
STD_IFG2STG_L_dp2ii = DCM.Ep.A{1,4}(3,5);
STD_IFG2STG_R_dp2ii = DCM.Ep.A{1,4}(4,6);

% MMN ones (conditional change): do not have a particular population
MMN_STG2A1_L = DCM.Ep.B{1,1}(1,3);
MMN_STG2A1_R = DCM.Ep.B{1,1}(2,4);
MMN_IFG2STG_L = DCM.Ep.B{1,1}(3,5);
MMN_IFG2STG_R = DCM.Ep.B{1,1}(4,6);

% Deep pyramidal to superficial pyramidal
DEV_STG2A1_L_dp2sp = DCM.Ep.A{1,3}(1,3) + MMN_STG2A1_L;
DEV_STG2A1_R_dp2sp = DCM.Ep.A{1,3}(2,4) + MMN_STG2A1_R;
DEV_IFG2STG_L_dp2sp = DCM.Ep.A{1,3}(3,5) + MMN_IFG2STG_L;
DEV_IFG2STG_R_dp2sp = DCM.Ep.A{1,3}(4,6) + MMN_IFG2STG_R;
% Deep pyramidal to inhibitory interneurons
DEV_STG2A1_L_dp2ii = DCM.Ep.A{1,4}(1,3) + MMN_STG2A1_L;
DEV_STG2A1_R_dp2ii = DCM.Ep.A{1,4}(2,4) + MMN_STG2A1_R;
DEV_IFG2STG_L_dp2ii = DCM.Ep.A{1,4}(3,5) + MMN_IFG2STG_L;
DEV_IFG2STG_R_dp2ii = DCM.Ep.A{1,4}(4,6) +  MMN_IFG2STG_R;

% Modulatory for the baseline/standard (M) and MMN (N): does not specify any neural populations
% MMN here refers to the conditional change in connections from baseline(standard) to deviant

STD_A1_L_mod = DCM.Ep.M(1,1);
STD_A1_R_mod = DCM.Ep.M(2,2);
STD_STG_L_mod = DCM.Ep.M(2,2);
STD_STG_R_mod = DCM.Ep.M(4,4);
STD_IFG_L_mod = DCM.Ep.M(5,5);
STD_IFG_R_mod = DCM.Ep.M(6,6);

MMN_A1_L_mod = DCM.Ep.N{1,1}(1,1);
MMN_A1_R_mod = DCM.Ep.N{1,1}(2,2);
MMN_STG_L_mod = DCM.Ep.N{1,1}(2,2);
MMN_STG_R_mod = DCM.Ep.N{1,1}(4,4);
MMN_IFG_L_mod = DCM.Ep.N{1,1}(5,5);
MMN_IFG_R_mod = DCM.Ep.N{1,1}(6,6);

DEV_A1_L_mod = STD_A1_L_mod + MMN_A1_L_mod;
DEV_A1_R_mod = STD_A1_R_mod + MMN_A1_R_mod;
DEV_STG_L_mod = STD_STG_L_mod + MMN_STG_L_mod;
DEV_STG_R_mod = STD_STG_R_mod + MMN_STG_R_mod;
DEV_IFG_L_mod = STD_IFG_L_mod + MMN_IFG_L_mod;
DEV_IFG_R_mod = STD_IFG_R_mod + MMN_IFG_R_mod;

% A and M objects: effective connectivity to the standard, 
% B and N objectS: conditional change from that standard in response to the deviant tones

%% Intrisinc connections within regions (CMC)

% Extract full from sparce matrix
G_matrix = full(DCM.Ep.G);

for roi = 1:length(ROIs)
    for ic = 1:length(intrinsic_connections)
        eval(['CMC_' ROIs{roi} '_' intrinsic_connections{ic} '= G_matrix(roi,ic);'])
    end
end

%% Add subject to structure for comparison
PEPP_MMN_Baseline_DCM{4,p+1} = STD_A1_L_ii;
PEPP_MMN_Baseline_DCM{5,p+1} = STD_A1_R_ii;
PEPP_MMN_Baseline_DCM{6,p+1} = STD_STG_L_ii;
PEPP_MMN_Baseline_DCM{7,p+1} = STD_STG_R_ii;
PEPP_MMN_Baseline_DCM{8,p+1} = STD_IFG_L_ii;
PEPP_MMN_Baseline_DCM{9,p+1} = STD_IFG_R_ii;
PEPP_MMN_Baseline_DCM{10,p+1} = STD_A1_L_ei;
PEPP_MMN_Baseline_DCM{11,p+1} = STD_A1_R_ei;
PEPP_MMN_Baseline_DCM{12,p+1} = STD_STG_L_ei;
PEPP_MMN_Baseline_DCM{13,p+1} = STD_STG_R_ei;
PEPP_MMN_Baseline_DCM{14,p+1} = STD_IFG_L_ei;
PEPP_MMN_Baseline_DCM{15,p+1} = STD_IFG_R_ei;
PEPP_MMN_Baseline_DCM{16,p+1} = STD_A1_L_py;
PEPP_MMN_Baseline_DCM{17,p+1} = STD_A1_R_py;
PEPP_MMN_Baseline_DCM{18,p+1} = STD_STG_L_py;
PEPP_MMN_Baseline_DCM{19,p+1} = STD_STG_R_py;
PEPP_MMN_Baseline_DCM{20,p+1} = STD_IFG_L_py;
PEPP_MMN_Baseline_DCM{21,p+1} = STD_IFG_R_py;
PEPP_MMN_Baseline_DCM{22,p+1} = DEV_A1_L_ii;
PEPP_MMN_Baseline_DCM{23,p+1} = DEV_A1_R_ii;
PEPP_MMN_Baseline_DCM{24,p+1} = DEV_STG_L_ii;
PEPP_MMN_Baseline_DCM{25,p+1} = DEV_STG_R_ii;
PEPP_MMN_Baseline_DCM{26,p+1} = DEV_IFG_L_ii;
PEPP_MMN_Baseline_DCM{27,p+1} = DEV_IFG_R_ii;
PEPP_MMN_Baseline_DCM{28,p+1} = DEV_A1_L_ei;
PEPP_MMN_Baseline_DCM{29,p+1} = DEV_A1_R_ei;
PEPP_MMN_Baseline_DCM{30,p+1} = DEV_STG_L_ei;
PEPP_MMN_Baseline_DCM{31,p+1} = DEV_STG_R_ei;
PEPP_MMN_Baseline_DCM{32,p+1} = DEV_IFG_L_ei;
PEPP_MMN_Baseline_DCM{33,p+1} = DEV_IFG_R_ei;
PEPP_MMN_Baseline_DCM{34,p+1} = DEV_A1_L_py;
PEPP_MMN_Baseline_DCM{35,p+1} = DEV_A1_R_py;
PEPP_MMN_Baseline_DCM{36,p+1} = DEV_STG_L_py;
PEPP_MMN_Baseline_DCM{37,p+1} = DEV_STG_R_py;
PEPP_MMN_Baseline_DCM{38,p+1} = DEV_IFG_L_py;
PEPP_MMN_Baseline_DCM{39,p+1} = DEV_IFG_R_py;

PEPP_MMN_Baseline_DCM{41,p+1} = STD_A12STG_L_sp2ss;
PEPP_MMN_Baseline_DCM{42,p+1} = STD_A12STG_R_sp2ss;
PEPP_MMN_Baseline_DCM{43,p+1} = STD_STG2IFG_L_sp2ss;
PEPP_MMN_Baseline_DCM{44,p+1} = STD_STG2IFG_R_sp2ss;
PEPP_MMN_Baseline_DCM{45,p+1} = STD_A12STG_L_sp2dp;
PEPP_MMN_Baseline_DCM{46,p+1} = STD_A12STG_R_sp2dp;
PEPP_MMN_Baseline_DCM{47,p+1} = STD_STG2IFG_L_sp2dp;
PEPP_MMN_Baseline_DCM{48,p+1} = STD_STG2IFG_R_sp2dp;
PEPP_MMN_Baseline_DCM{49,p+1} = DEV_A12STG_L_sp2ss;
PEPP_MMN_Baseline_DCM{50,p+1} = DEV_A12STG_R_sp2ss;
PEPP_MMN_Baseline_DCM{51,p+1} = DEV_STG2IFG_L_sp2ss;
PEPP_MMN_Baseline_DCM{52,p+1} = DEV_STG2IFG_R_sp2ss;
PEPP_MMN_Baseline_DCM{53,p+1} = DEV_A12STG_L_sp2dp;
PEPP_MMN_Baseline_DCM{54,p+1} = DEV_A12STG_R_sp2dp;
PEPP_MMN_Baseline_DCM{55,p+1} = DEV_STG2IFG_L_sp2dp;
PEPP_MMN_Baseline_DCM{56,p+1} = DEV_STG2IFG_R_sp2dp;
PEPP_MMN_Baseline_DCM{57,p+1} = MMN_A12STG_L;
PEPP_MMN_Baseline_DCM{58,p+1} = MMN_A12STG_R;
PEPP_MMN_Baseline_DCM{59,p+1} = MMN_STG2IFG_L;
PEPP_MMN_Baseline_DCM{60,p+1} = MMN_STG2IFG_R;

PEPP_MMN_Baseline_DCM{62,p+1} = STD_STG2A1_L_dp2sp;
PEPP_MMN_Baseline_DCM{63,p+1} = STD_STG2A1_R_dp2sp;
PEPP_MMN_Baseline_DCM{64,p+1} = STD_IFG2STG_L_dp2sp;
PEPP_MMN_Baseline_DCM{65,p+1} = STD_IFG2STG_R_dp2sp;
PEPP_MMN_Baseline_DCM{66,p+1} = STD_STG2A1_L_dp2ii;
PEPP_MMN_Baseline_DCM{67,p+1} = STD_STG2A1_R_dp2ii;
PEPP_MMN_Baseline_DCM{68,p+1} = STD_IFG2STG_L_dp2ii;
PEPP_MMN_Baseline_DCM{69,p+1} = STD_IFG2STG_R_dp2ii;
PEPP_MMN_Baseline_DCM{70,p+1} = DEV_STG2A1_L_dp2sp;
PEPP_MMN_Baseline_DCM{71,p+1} = DEV_STG2A1_R_dp2sp;
PEPP_MMN_Baseline_DCM{72,p+1} = DEV_IFG2STG_L_dp2sp;
PEPP_MMN_Baseline_DCM{73,p+1} = DEV_IFG2STG_R_dp2sp;
PEPP_MMN_Baseline_DCM{74,p+1} = DEV_STG2A1_L_dp2ii;
PEPP_MMN_Baseline_DCM{75,p+1} = DEV_STG2A1_R_dp2ii;
PEPP_MMN_Baseline_DCM{76,p+1} = DEV_IFG2STG_L_dp2ii;
PEPP_MMN_Baseline_DCM{77,p+1} = DEV_IFG2STG_R_dp2ii;
PEPP_MMN_Baseline_DCM{78,p+1} = MMN_STG2A1_L;
PEPP_MMN_Baseline_DCM{79,p+1} = MMN_STG2A1_R;
PEPP_MMN_Baseline_DCM{80,p+1} = MMN_IFG2STG_L;
PEPP_MMN_Baseline_DCM{81,p+1} = MMN_IFG2STG_R;

PEPP_MMN_Baseline_DCM{83,p+1} = STD_A1_L_mod;
PEPP_MMN_Baseline_DCM{84,p+1} = STD_A1_R_mod;
PEPP_MMN_Baseline_DCM{85,p+1} = STD_STG_L_mod;
PEPP_MMN_Baseline_DCM{86,p+1} = STD_STG_R_mod;
PEPP_MMN_Baseline_DCM{87,p+1} = STD_IFG_L_mod;
PEPP_MMN_Baseline_DCM{88,p+1} = STD_IFG_R_mod;
PEPP_MMN_Baseline_DCM{89,p+1} = DEV_A1_L_mod;
PEPP_MMN_Baseline_DCM{90,p+1} = DEV_A1_R_mod;
PEPP_MMN_Baseline_DCM{91,p+1} = DEV_STG_L_mod;
PEPP_MMN_Baseline_DCM{92,p+1} = DEV_STG_R_mod;
PEPP_MMN_Baseline_DCM{93,p+1} = DEV_IFG_L_mod;
PEPP_MMN_Baseline_DCM{94,p+1} = DEV_IFG_R_mod;
PEPP_MMN_Baseline_DCM{95,p+1} = MMN_A1_L_mod;
PEPP_MMN_Baseline_DCM{96,p+1} = MMN_A1_R_mod;
PEPP_MMN_Baseline_DCM{97,p+1} = MMN_STG_L_mod;
PEPP_MMN_Baseline_DCM{98,p+1} = MMN_STG_R_mod;
PEPP_MMN_Baseline_DCM{99,p+1} = MMN_IFG_L_mod;
PEPP_MMN_Baseline_DCM{100,p+1} = MMN_IFG_R_mod;

PEPP_MMN_Baseline_DCM{102,p+1} = CMC_A1_L_ii2ss;
PEPP_MMN_Baseline_DCM{103,p+1} = CMC_A1_L_sp2sp;
PEPP_MMN_Baseline_DCM{104,p+1} = CMC_A1_L_sp2ss;
PEPP_MMN_Baseline_DCM{105,p+1} = CMC_A1_R_ii2ss;
PEPP_MMN_Baseline_DCM{106,p+1} = CMC_A1_R_sp2sp;
PEPP_MMN_Baseline_DCM{107,p+1} = CMC_A1_R_sp2ss;
PEPP_MMN_Baseline_DCM{108,p+1} = CMC_STG_L_ii2ss;
PEPP_MMN_Baseline_DCM{109,p+1} = CMC_STG_L_sp2sp;
PEPP_MMN_Baseline_DCM{110,p+1} = CMC_STG_L_sp2ss;
PEPP_MMN_Baseline_DCM{111,p+1} = CMC_STG_R_ii2ss;
PEPP_MMN_Baseline_DCM{112,p+1} = CMC_STG_R_sp2sp;
PEPP_MMN_Baseline_DCM{113,p+1} = CMC_STG_R_sp2ss;
PEPP_MMN_Baseline_DCM{114,p+1} = CMC_IFG_L_ii2ss;
PEPP_MMN_Baseline_DCM{115,p+1} = CMC_IFG_L_sp2sp;
PEPP_MMN_Baseline_DCM{116,p+1} = CMC_IFG_L_sp2ss;
PEPP_MMN_Baseline_DCM{117,p+1} = CMC_IFG_R_ii2ss;
PEPP_MMN_Baseline_DCM{118,p+1} = CMC_IFG_R_sp2sp;
PEPP_MMN_Baseline_DCM{119,p+1} = CMC_IFG_R_sp2ss;

end

%% Add average and std dev/err of measures (separately for C and FE) and save

last_column = size(PEPP_MMN_Baseline_DCM,2) + 1;
PEPP_MMN_Baseline_DCM{1,last_column} = 'AVERAGE_HC'; column_average_c = last_column;
PEPP_MMN_Baseline_DCM{1,last_column + 1} = 'AVERAGE_FE'; column_average_fe = last_column + 1;
PEPP_MMN_Baseline_DCM{1,last_column + 2} = 'STD_DEV_HC'; column_std_dev_c = last_column + 2;
PEPP_MMN_Baseline_DCM{1,last_column + 3} = 'STD_DEV_FE'; column_std_dev_fe = last_column + 3;
PEPP_MMN_Baseline_DCM{1,last_column + 4} = 'STD_ERR_HC'; column_std_err_c = last_column + 4;
PEPP_MMN_Baseline_DCM{1,last_column + 5} = 'STD_ERR_FE'; column_std_err_fe = last_column + 5;

% Source waveforms first
rows_to_average = 4:39;
for ca = 1:length(rows_to_average) 
PEPP_MMN_Baseline_DCM{rows_to_average(ca),column_average_fe} = mean(cell2mat(PEPP_MMN_Baseline_DCM(rows_to_average(ca),pos_fe)),2);
PEPP_MMN_Baseline_DCM{rows_to_average(ca),column_average_c} = mean(cell2mat(PEPP_MMN_Baseline_DCM(rows_to_average(ca),pos_c)),2);
PEPP_MMN_Baseline_DCM{rows_to_average(ca),column_std_dev_fe} = std(cell2mat(PEPP_MMN_Baseline_DCM(rows_to_average(ca),pos_fe)),0,2);
PEPP_MMN_Baseline_DCM{rows_to_average(ca),column_std_dev_c} = std(cell2mat(PEPP_MMN_Baseline_DCM(rows_to_average(ca),pos_c)),0,2);
PEPP_MMN_Baseline_DCM{rows_to_average(ca),column_std_err_fe} = PEPP_MMN_Baseline_DCM{rows_to_average(ca),column_std_dev_fe}/sqrt(length(pos_fe));
PEPP_MMN_Baseline_DCM{rows_to_average(ca),column_std_err_c} = PEPP_MMN_Baseline_DCM{rows_to_average(ca),column_std_dev_c}/sqrt(length(pos_c));
end

% Rest of measures
rows_to_average = [41:60 62:81 83:100 102:119];
for ca = 1:length(rows_to_average)
PEPP_MMN_Baseline_DCM{rows_to_average(ca),column_average_fe} = mean(cell2mat(PEPP_MMN_Baseline_DCM(rows_to_average(ca),pos_fe)));
PEPP_MMN_Baseline_DCM{rows_to_average(ca),column_average_c} = mean(cell2mat(PEPP_MMN_Baseline_DCM(rows_to_average(ca),pos_c)));
end

% Save pos_fe and pos_c for when we plot using these
save([outcome_path '/Results_DCM/pos_subj_CMC_' MMN_types{tymmn} '_' modality_data{mode} '.mat'],'pos_c','pos_fe');

% Save matrix 
save([outcome_path '/Results_DCM/DCM_CMC_' MMN_types{tymmn} '_' modality_data{mode} '.mat'],'PEPP_MMN_Baseline_DCM');

disp(['Done extracting DCM CMC values (' MMN_types{tymmn} '_' modality_data{mode} ')']);

end

end

end

%% Define Plot-specific variables

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SPECIFIC PARAMETERS OF THIS SECTION %%%%%%%%%%%%%%%%%%%%%%
set(0,'defaultfigurecolor',[1 1 1]); % I want white backgrounds
participant_group = {'FE','C'};
populations_waveform = {'ii','ei','py'};
ROIs2plot = {'A1','STG','IFG'};
hemisphere = {'L','R'};
dev_GAVR = 2; % 1 = Standard Deviation 2 = Standard Error
plot_baseline = -50;
plot_post = 300;
s_r = 1000;
time_samples=linspace(plot_baseline,plot_post,((((plot_baseline*(-1)) + plot_post)/1000)*s_r) +1);
% Colors
color_group = {[0 220 0],[0 0 0]}; % FE C
color_hem = {[0 0 255],[255 0 0]}; % L R
color_roi = {[255 0 255],[255 0 0],[150 0 0],...
    [255 128 0],[0 153 76],[0 0 153]}; % Depends on how many you put above (bot only one hem)
transparency = 0.2;
shaded_areas = 0; % 0 = NO shaded areas; 1 = YES
shaded_areas_scalp_Pitch_MMN_EEG = {[70 170 170 70]}; % {[27 33 33 27],[51 61 61 51]};;
shaded_areas_scalp_Pitch_MMN_MEG = {[70 170 170 70]}; % {[27 33 33 27],[51 61 61 51]};;
shaded_areas_scalp_Dur_MMN_EEG = {[110 210 210 110]};
shaded_areas_scalp_Dur_MMN_MEG = {[110 210 210 110]};
shaded_areas_source_left_Pitch_MMN = {[70 170 170 70]}; % {[110 130 130 110]};
shaded_areas_source_left_Dur_MMN = {[110 210 210 110]};
shaded_areas_source_right_Pitch_MMN = {[70 170 170 70]}; % {[27 33 33 27],[51 61 61 51]};
shaded_areas_source_right_Dur_MMN = {[110 210 210 110]};
save_figures_option = 0; % 0 = NO 1 = YES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot butterfly source waveforms DCM 

if strcmp(plot_butterfly,'YES')
for mode = 1:length(modality_data)
for tymmn = 1:length(MMN_types)
    
load([outcome_path '/Results_DCM/DCM_CMC_' MMN_types{tymmn} '_' modality_data{mode} '.mat']);
load([outcome_path '/Results_DCM/pos_subj_CMC_' MMN_types{tymmn} '_' modality_data{mode} '.mat']);

cond_titles = {'STD','DEV'};
if strcmp(MMN_types{tymmn},'pMMN')
    title_plot = {'STD','PDEV'};
elseif strcmp(MMN_types{tymmn},'dMMN')
    title_plot = {'STD','DDEV'};
end

for ct = 1:length(cond_titles)
for pop = 1:length(populations_waveform)
    if strcmp(populations_waveform{pop},'ii')
        line_style = ':';
    elseif strcmp(populations_waveform{pop},'ei')
        line_style = '-.';
    elseif strcmp(populations_waveform{pop},'py')
        line_style = '-';
    end
    % New plot 
    figure;
    % figure('units','normalized','outerposition',[0 0 1 1]);
    num_subplots = 12; % For later
    for i = 1:num_subplots
        h(i) = subplot(3,4,i);
    end
for roi = 1:length(ROIs2plot)
    current_max = []; current_min = []; pos_ylim = 1; % To adjust ylim later
for hem = 1:length(hemisphere)    
for pg = 1:length(participant_group)
    if strcmp(participant_group{pg},'FE')
        participant_list = PEPP_MMN_Baseline_DCM(1,pos_fe);
    elseif strcmp(participant_group{pg},'C')
        participant_list = PEPP_MMN_Baseline_DCM(1,pos_c);
    end

    % Determine position in figure based on iteration
    if strcmp(ROIs2plot{roi},'A1') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'C') 
        plot_pos = 1;
    elseif strcmp(ROIs2plot{roi},'A1') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'C') 
        plot_pos = 2;
    elseif strcmp(ROIs2plot{roi},'A1') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'FE') 
        plot_pos = 3;
    elseif strcmp(ROIs2plot{roi},'A1') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'FE') 
        plot_pos = 4;
    elseif strcmp(ROIs2plot{roi},'STG') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'C') 
        plot_pos = 5;
    elseif strcmp(ROIs2plot{roi},'STG') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'C') 
        plot_pos = 6;
    elseif strcmp(ROIs2plot{roi},'STG') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'FE') 
        plot_pos = 7;
    elseif strcmp(ROIs2plot{roi},'STG') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'FE') 
        plot_pos = 8;
    elseif strcmp(ROIs2plot{roi},'IFG') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'C') 
        plot_pos = 9;
    elseif strcmp(ROIs2plot{roi},'IFG') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'C') 
        plot_pos = 10;
    elseif strcmp(ROIs2plot{roi},'IFG') && strcmp(hemisphere{hem},'L') && strcmp(participant_group{pg},'FE') 
        plot_pos = 11;
    elseif strcmp(ROIs2plot{roi},'IFG') && strcmp(hemisphere{hem},'R') && strcmp(participant_group{pg},'FE') 
        plot_pos = 12;
    end
    row_to_look = find(strcmp(PEPP_MMN_Baseline_DCM(:,1),[cond_titles{ct} '_' ROIs2plot{roi} '_' hemisphere{hem} '_' populations_waveform{pop}]));
    
    % Now plot
    for p = 1:length(participant_list)
        if strcmp(participant_group{pg},'FE')
            pos_subject = pos_fe(p);
         elseif strcmp(participant_group{pg},'C')
            pos_subject = pos_c(p);
         end
        % If cluster was selected, average amplitudes before plotting
        Amplitudes = PEPP_MMN_Baseline_DCM(row_to_look,pos_subject);
        hold (h(plot_pos),'on')
        plot(h(plot_pos),time_samples, Amplitudes{:},line_style);
        current_max(pos_ylim) = max(Amplitudes{:});
        current_min(pos_ylim) = min(Amplitudes{:});
        pos_ylim = pos_ylim + 1; % iterations of participants and conditions
    end
    % Add title (unique of every subplot) and legend (unique of group)
    hold (h(plot_pos),'on')
    title (h(plot_pos),[ROIs2plot{roi} ' ' hemisphere{hem} ' ' participant_group{pg}]);
    legend(h(plot_pos),participant_list);
end
end
    if roi == 1
        which_subplots = 1:4;
    elseif roi == 2
        which_subplots = 5:8;
    elseif roi == 3
        which_subplots = 9:12;
    end
    % Once data is put in figure, adjust format
    y_lim_min = min(current_min);
    y_lim_max = max(current_max);
    for i = which_subplots
        hold (h(i),'on')
        xlim(h(i),[plot_baseline,plot_post]);
        ylim(h(i),[y_lim_min,y_lim_max]);
        plot(h(i),xlim,[0 0], '-k')
        line(h(i),[0 0], [y_lim_min y_lim_max],'Color','black')
    end
end 

    hold (h(9),'on')
    ylabel(h(9),'Relative units'); 
    xlabel(h(9),'Time (ms)');
    current_title = [modality_data{mode} '_' title_plot{ct} '_' populations_waveform{pop}];
    current_title = strrep(current_title,'_',' ');
    suptitle(current_title)
    pom=findobj('type','legend');
    delete(pom);

    if shaded_areas == 1
        if strcmp(modality_data{mode},'EEG')
            shaded_areas_scalp_Pitch_MMN = shaded_areas_scalp_Pitch_MMN_EEG;
            shaded_areas_scalp_Dur_MMN = shaded_areas_scalp_Dur_MMN_EEG;
        elseif strcmp(modality_data{mode},'MEG')
            shaded_areas_scalp_Pitch_MMN = shaded_areas_scalp_Pitch_MMN_MEG;
            shaded_areas_scalp_Dur_MMN = shaded_areas_scalp_Dur_MMN_MEG; 
        end
    % Patch Pitch MMN plots
    for i = [4 9] 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_scalp_Pitch_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_scalp_Pitch_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    % Patch Duration MMN plots
    for i = [5 10] 
        hold (h(i),'on')
        for shad = 1:length(shaded_areas_scalp_Dur_MMN)
            gray = [0 0 0];
            patch(h(i),shaded_areas_scalp_Dur_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
        end
    end
    end

end
end

end
end
end

%% Plot GAVR sources DCM

if strcmp(plot_gavr,'YES')
for mode = 1:length(modality_data)
for tymmn = 1:length(MMN_types)
    
load([outcome_path '/Results_DCM/DCM_CMC_' MMN_types{tymmn} '_' modality_data{mode} '.mat']);
load([outcome_path '/Results_DCM/pos_subj_CMC_' MMN_types{tymmn} '_' modality_data{mode} '.mat']);

cond_titles = {'STD','DEV'};
if strcmp(MMN_types{tymmn},'pMMN')
    title_plot = {'STD','PDEV'};
elseif strcmp(MMN_types{tymmn},'dMMN')
    title_plot = {'STD','DDEV'};
end
    
% New plot 
figure;
% figure('units','normalized','outerposition',[0 0 1 1]);
num_subplots = 36; % For later
for i = 1:num_subplots
    h(i) = subplot(6,6,i);
end

for roi = 1:length(ROIs)
for pop = 1:length(populations_waveform)
    if strcmp(populations_waveform{pop},'ii')
        line_style = ':';
    elseif strcmp(populations_waveform{pop},'ei')
        line_style = '-.';
    elseif strcmp(populations_waveform{pop},'py')
        line_style = '-';
    end
    current_max = []; current_min = []; pos_ylim = 1; % To adjust ylim later
for ct = 1:length(cond_titles)  
% Determine position in figure based on iteration
if strcmp(ROIs{roi},'A1_L') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'STD') 
    plot_pos = 1;
elseif strcmp(ROIs{roi},'A1_L') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 2;
elseif strcmp(ROIs{roi},'A1_L') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'STD')
    plot_pos = 3;
elseif strcmp(ROIs{roi},'A1_L') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 4;
elseif strcmp(ROIs{roi},'A1_L') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'STD')
    plot_pos = 5;
elseif strcmp(ROIs{roi},'A1_L') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 6;
elseif strcmp(ROIs{roi},'A1_R') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'STD') 
    plot_pos = 7;
elseif strcmp(ROIs{roi},'A1_R') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 8;
elseif strcmp(ROIs{roi},'A1_R') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'STD')
    plot_pos = 9;
elseif strcmp(ROIs{roi},'A1_R') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 10;
elseif strcmp(ROIs{roi},'A1_R') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'STD')
    plot_pos = 11;
elseif strcmp(ROIs{roi},'A1_R') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 12;
elseif strcmp(ROIs{roi},'STG_L') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'STD') 
    plot_pos = 13;
elseif strcmp(ROIs{roi},'STG_L') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 14;
elseif strcmp(ROIs{roi},'STG_L') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'STD')
    plot_pos = 15;
elseif strcmp(ROIs{roi},'STG_L') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 16;
elseif strcmp(ROIs{roi},'STG_L') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'STD')
    plot_pos = 17;
elseif strcmp(ROIs{roi},'STG_L') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 18;
elseif strcmp(ROIs{roi},'STG_R') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'STD') 
    plot_pos = 19;
elseif strcmp(ROIs{roi},'STG_R') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 20;
elseif strcmp(ROIs{roi},'STG_R') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'STD')
    plot_pos = 21;
elseif strcmp(ROIs{roi},'STG_R') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 22;
elseif strcmp(ROIs{roi},'STG_R') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'STD')
    plot_pos = 23;
elseif strcmp(ROIs{roi},'STG_R') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 24;
elseif strcmp(ROIs{roi},'IFG_L') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'STD') 
    plot_pos = 25;
elseif strcmp(ROIs{roi},'IFG_L') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 26;
elseif strcmp(ROIs{roi},'IFG_L') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'STD')
    plot_pos = 27;
elseif strcmp(ROIs{roi},'IFG_L') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 28;
elseif strcmp(ROIs{roi},'IFG_L') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'STD')
    plot_pos = 29;
elseif strcmp(ROIs{roi},'IFG_L') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 30;
elseif strcmp(ROIs{roi},'IFG_R') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'STD') 
    plot_pos = 31;
elseif strcmp(ROIs{roi},'IFG_R') && strcmp(populations_waveform{pop},'ii') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 32;
elseif strcmp(ROIs{roi},'IFG_R') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'STD')
    plot_pos = 33;
elseif strcmp(ROIs{roi},'IFG_R') && strcmp(populations_waveform{pop},'ei') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 34;
elseif strcmp(ROIs{roi},'IFG_R') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'STD')
    plot_pos = 35;
elseif strcmp(ROIs{roi},'IFG_R') && strcmp(populations_waveform{pop},'py') && strcmp(cond_titles{ct},'DEV')
    plot_pos = 36;
end
row_to_look = find(strcmp(PEPP_MMN_Baseline_DCM(:,1),[cond_titles{ct} '_' ROIs{roi} '_' populations_waveform{pop}]));

% Now plot
for pg = 1:length(participant_group)
    column_to_look = strcmp(PEPP_MMN_Baseline_DCM(1,:),'AVERAGE_HC');
    if strcmp(participant_group{pg},'FE')
        column_average = strcmp(PEPP_MMN_Baseline_DCM(1,:),'AVERAGE_FE');
        column_std_dev = strcmp(PEPP_MMN_Baseline_DCM(1,:),'STD_DEV_FE');
        column_std_err = strcmp(PEPP_MMN_Baseline_DCM(1,:),'STD_ERR_FE');
     elseif strcmp(participant_group{pg},'C')
        column_average = strcmp(PEPP_MMN_Baseline_DCM(1,:),'AVERAGE_HC');
        column_std_dev = strcmp(PEPP_MMN_Baseline_DCM(1,:),'STD_DEV_HC');
        column_std_err = strcmp(PEPP_MMN_Baseline_DCM(1,:),'STD_ERR_HC');
     end
    % If cluster was selected, average amplitudes before plotting
    average = PEPP_MMN_Baseline_DCM(row_to_look,column_average); 
    average = average{:}';
    
    % Load stdev or stderr
    if dev_GAVR == 1 % standard deviation
        dev = PEPP_MMN_Baseline_DCM(row_to_look,column_std_dev);
    elseif dev_GAVR == 2 % standard error
        dev = PEPP_MMN_Baseline_DCM(row_to_look,column_std_err);
    end  
    dev = dev{:}';
    
    % Set data ready for plot
    curve1 = average + dev;
    curve2 = average - dev;
    time_samples_2 = [time_samples, fliplr(time_samples)];
    inBetween = [curve1, fliplr(curve2)];
    % Also, grab values for ylim later
    current_max(pos_ylim) = max(curve1);
    current_min(pos_ylim) = min(curve2);
    pos_ylim = pos_ylim + 1; % iterations of participants and conditions
    
    % Now plot
    hold (h(plot_pos),'on')
    fill(h(plot_pos),time_samples_2, inBetween, (color_group{pg}/256), 'FaceAlpha', transparency, 'LineStyle', 'none','HandleVisibility','off');
    plot(h(plot_pos),time_samples, average, 'color', (color_group{pg}/256), 'LineWidth', 1.5, 'linestyle', line_style);
end
% Add title (unique of every subplot) and legend (unique of group)
hold (h(plot_pos),'on')
this_title = [ROIs{roi} ' ' populations_waveform{pop} ' ' cond_titles{ct}];
this_title = strrep(this_title,'_',' ');
title(h(plot_pos),this_title);
end 
if roi == 1 && pop == 1
    which_subplots = 1:2;
elseif roi == 1 && pop == 2
    which_subplots = 3:4;
elseif roi == 1 && pop == 3
    which_subplots = 5:6;
elseif roi == 2 && pop == 1
    which_subplots = 7:8;
elseif roi == 2 && pop == 2
    which_subplots = 9:10;
elseif roi == 2 && pop == 3
    which_subplots = 11:12;
elseif roi == 3 && pop == 1
    which_subplots = 13:14;
elseif roi == 3 && pop == 2
    which_subplots = 15:16;
elseif roi == 3 && pop == 3
    which_subplots = 17:18;
elseif roi == 4 && pop == 1
    which_subplots = 19:20;
elseif roi == 4 && pop == 2
    which_subplots = 21:22;
elseif roi == 4 && pop == 3
    which_subplots = 23:24;
elseif roi == 5 && pop == 1
    which_subplots = 25:26;
elseif roi == 5 && pop == 2
    which_subplots = 27:28;
elseif roi == 5 && pop == 3
    which_subplots = 29:30;
elseif roi == 6 && pop == 1
    which_subplots = 31:32;
elseif roi == 6 && pop == 2
    which_subplots = 33:34;
elseif roi == 6 && pop == 3
    which_subplots = 35:36;
end

% Once data is put in figure, adjust format
y_lim_min = min(current_min);
y_lim_max = max(current_max);
for i = which_subplots
    hold (h(i),'on')
    ylim(h(i),[y_lim_min,y_lim_max]);
    plot(h(i),xlim,[0 0], '-k')
    line(h(i),[0 0], [y_lim_min y_lim_max],'Color','black')
end
end
end

hold (h(1),'on')
legend(h(1),participant_group);
hold (h(31),'on')
ylabel(h(31),'Relative units'); 
xlabel(h(31),'Time (ms)');
for i = 1:num_subplots
    hold (h(i),'on')
    xlim(h(i),[plot_baseline,plot_post]);
end
for i = 1:31
    hold (h(i),'on')
    set(h(i),'XTickLabel',[]);
end
current_title = [modality_data{mode} ' source solutions DCM (' MMN_types{tymmn} ')'];
current_title = strrep(current_title,'_',' ');
suptitle(current_title)

if shaded_areas == 1
    if strcmp(modality_data{mode},'EEG')
        shaded_areas_scalp_Pitch_MMN = shaded_areas_scalp_Pitch_MMN_EEG;
        shaded_areas_scalp_Dur_MMN = shaded_areas_scalp_Dur_MMN_EEG;
    elseif strcmp(modality_data{mode},'MEG')
        shaded_areas_scalp_Pitch_MMN = shaded_areas_scalp_Pitch_MMN_MEG;
        shaded_areas_scalp_Dur_MMN = shaded_areas_scalp_Dur_MMN_MEG; 
    end
% Patch Pitch MMN plots
for i = [4 9] 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_scalp_Pitch_MMN)
        gray = [0 0 0];
        patch(h(i),shaded_areas_scalp_Pitch_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
% Patch Duration MMN plots
for i = [5 10] 
    hold (h(i),'on')
    for shad = 1:length(shaded_areas_scalp_Dur_MMN)
        gray = [0 0 0];
        patch(h(i),shaded_areas_scalp_Dur_MMN{shad},[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', transparency,'EdgeAlpha',0.1,'HandleVisibility','off')
    end
end
end

end
end
end

%% Plot connectivity matrices

if strcmp(plot_connectivity,'YES')
for mode = 1:length(modality_data)
for tymmn = 1:length(MMN_types)
    
load([outcome_path '/Results_DCM/DCM_CMC_' MMN_types{tymmn} '_' modality_data{mode} '.mat']);
load([outcome_path '/Results_DCM/pos_subj_CMC_' MMN_types{tymmn} '_' modality_data{mode} '.mat']);

cond_titles = {'STD','DEV'};
if strcmp(MMN_types{tymmn},'pMMN')
    title_plot = {'STD','PDEV'};
elseif strcmp(MMN_types{tymmn},'dMMN')
    title_plot = {'STD','DDEV'};
end

avg_c_pos = find(strcmp(PEPP_MMN_Baseline_DCM(1,:),'AVERAGE_HC'));
avg_fe_pos = find(strcmp(PEPP_MMN_Baseline_DCM(1,:),'AVERAGE_FE'));

% Define connectivity matrices
if 1 == 1 % silly way to compress it
STD_sp2ss_dp2sp_matrix_C = [];
DEV_sp2ss_dp2sp_matrix_C = [];
STD_sp2dp_dp2ii_matrix_C = [];
DEV_sp2dp_dp2ii_matrix_C = [];
MMN_matrix_C = [];

STD_sp2ss_dp2sp_matrix_FE = [];
DEV_sp2ss_dp2sp_matrix_FE = [];
STD_sp2dp_dp2ii_matrix_FE = [];
DEV_sp2dp_dp2ii_matrix_FE = [];
MMN_matrix_FE = [];

STD_sp2ss_dp2sp_matrix_C(3,1) = cell2mat(PEPP_MMN_Baseline_DCM(41,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(4,2) = cell2mat(PEPP_MMN_Baseline_DCM(42,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(5,3) = cell2mat(PEPP_MMN_Baseline_DCM(43,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(6,4) = cell2mat(PEPP_MMN_Baseline_DCM(44,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(1,3) = cell2mat(PEPP_MMN_Baseline_DCM(62,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(2,4) = cell2mat(PEPP_MMN_Baseline_DCM(63,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(3,5) = cell2mat(PEPP_MMN_Baseline_DCM(64,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(4,6) = cell2mat(PEPP_MMN_Baseline_DCM(65,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(1,1) = cell2mat(PEPP_MMN_Baseline_DCM(83,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(2,2) = cell2mat(PEPP_MMN_Baseline_DCM(84,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(3,3) = cell2mat(PEPP_MMN_Baseline_DCM(85,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(4,4) = cell2mat(PEPP_MMN_Baseline_DCM(86,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(5,5) = cell2mat(PEPP_MMN_Baseline_DCM(87,avg_c_pos));
STD_sp2ss_dp2sp_matrix_C(6,6) = cell2mat(PEPP_MMN_Baseline_DCM(88,avg_c_pos));

STD_sp2dp_dp2ii_matrix_C(3,1) = cell2mat(PEPP_MMN_Baseline_DCM(45,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(4,2) = cell2mat(PEPP_MMN_Baseline_DCM(46,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(5,3) = cell2mat(PEPP_MMN_Baseline_DCM(47,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(6,4) = cell2mat(PEPP_MMN_Baseline_DCM(48,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(1,3) = cell2mat(PEPP_MMN_Baseline_DCM(66,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(2,4) = cell2mat(PEPP_MMN_Baseline_DCM(67,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(3,5) = cell2mat(PEPP_MMN_Baseline_DCM(68,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(4,6) = cell2mat(PEPP_MMN_Baseline_DCM(69,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(1,1) = cell2mat(PEPP_MMN_Baseline_DCM(83,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(2,2) = cell2mat(PEPP_MMN_Baseline_DCM(84,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(3,3) = cell2mat(PEPP_MMN_Baseline_DCM(85,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(4,4) = cell2mat(PEPP_MMN_Baseline_DCM(86,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(5,5) = cell2mat(PEPP_MMN_Baseline_DCM(87,avg_c_pos));
STD_sp2dp_dp2ii_matrix_C(6,6) = cell2mat(PEPP_MMN_Baseline_DCM(88,avg_c_pos));

DEV_sp2ss_dp2sp_matrix_C(3,1) = cell2mat(PEPP_MMN_Baseline_DCM(49,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(4,2) = cell2mat(PEPP_MMN_Baseline_DCM(50,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(5,3) = cell2mat(PEPP_MMN_Baseline_DCM(51,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(6,4) = cell2mat(PEPP_MMN_Baseline_DCM(52,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(1,3) = cell2mat(PEPP_MMN_Baseline_DCM(70,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(2,4) = cell2mat(PEPP_MMN_Baseline_DCM(71,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(3,5) = cell2mat(PEPP_MMN_Baseline_DCM(72,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(4,6) = cell2mat(PEPP_MMN_Baseline_DCM(73,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(1,1) = cell2mat(PEPP_MMN_Baseline_DCM(89,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(2,2) = cell2mat(PEPP_MMN_Baseline_DCM(90,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(3,3) = cell2mat(PEPP_MMN_Baseline_DCM(91,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(4,4) = cell2mat(PEPP_MMN_Baseline_DCM(92,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(5,5) = cell2mat(PEPP_MMN_Baseline_DCM(93,avg_c_pos));
DEV_sp2ss_dp2sp_matrix_C(6,6) = cell2mat(PEPP_MMN_Baseline_DCM(94,avg_c_pos));

DEV_sp2dp_dp2ii_matrix_C(3,1) = cell2mat(PEPP_MMN_Baseline_DCM(53,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(4,2) = cell2mat(PEPP_MMN_Baseline_DCM(54,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(5,3) = cell2mat(PEPP_MMN_Baseline_DCM(55,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(6,4) = cell2mat(PEPP_MMN_Baseline_DCM(56,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(1,3) = cell2mat(PEPP_MMN_Baseline_DCM(74,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(2,4) = cell2mat(PEPP_MMN_Baseline_DCM(75,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(3,5) = cell2mat(PEPP_MMN_Baseline_DCM(76,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(4,6) = cell2mat(PEPP_MMN_Baseline_DCM(77,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(1,1) = cell2mat(PEPP_MMN_Baseline_DCM(89,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(2,2) = cell2mat(PEPP_MMN_Baseline_DCM(90,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(3,3) = cell2mat(PEPP_MMN_Baseline_DCM(91,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(4,4) = cell2mat(PEPP_MMN_Baseline_DCM(92,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(5,5) = cell2mat(PEPP_MMN_Baseline_DCM(93,avg_c_pos));
DEV_sp2dp_dp2ii_matrix_C(6,6) = cell2mat(PEPP_MMN_Baseline_DCM(94,avg_c_pos));

MMN_matrix_C(3,1) = cell2mat(PEPP_MMN_Baseline_DCM(57,avg_c_pos));
MMN_matrix_C(4,2) = cell2mat(PEPP_MMN_Baseline_DCM(58,avg_c_pos));
MMN_matrix_C(5,3) = cell2mat(PEPP_MMN_Baseline_DCM(59,avg_c_pos));
MMN_matrix_C(6,4) = cell2mat(PEPP_MMN_Baseline_DCM(60,avg_c_pos));
MMN_matrix_C(1,3) = cell2mat(PEPP_MMN_Baseline_DCM(78,avg_c_pos));
MMN_matrix_C(2,4) = cell2mat(PEPP_MMN_Baseline_DCM(79,avg_c_pos));
MMN_matrix_C(3,5) = cell2mat(PEPP_MMN_Baseline_DCM(80,avg_c_pos));
MMN_matrix_C(4,6) = cell2mat(PEPP_MMN_Baseline_DCM(81,avg_c_pos));
MMN_matrix_C(1,1) = cell2mat(PEPP_MMN_Baseline_DCM(95,avg_c_pos));
MMN_matrix_C(2,2) = cell2mat(PEPP_MMN_Baseline_DCM(96,avg_c_pos));
MMN_matrix_C(3,3) = cell2mat(PEPP_MMN_Baseline_DCM(97,avg_c_pos));
MMN_matrix_C(4,4) = cell2mat(PEPP_MMN_Baseline_DCM(98,avg_c_pos));
MMN_matrix_C(5,5) = cell2mat(PEPP_MMN_Baseline_DCM(99,avg_c_pos));
MMN_matrix_C(6,6) = cell2mat(PEPP_MMN_Baseline_DCM(100,avg_c_pos));

STD_sp2ss_dp2sp_matrix_FE(3,1) = cell2mat(PEPP_MMN_Baseline_DCM(41,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(4,2) = cell2mat(PEPP_MMN_Baseline_DCM(42,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(5,3) = cell2mat(PEPP_MMN_Baseline_DCM(43,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(6,4) = cell2mat(PEPP_MMN_Baseline_DCM(44,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(1,3) = cell2mat(PEPP_MMN_Baseline_DCM(62,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(2,4) = cell2mat(PEPP_MMN_Baseline_DCM(63,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(3,5) = cell2mat(PEPP_MMN_Baseline_DCM(64,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(4,6) = cell2mat(PEPP_MMN_Baseline_DCM(65,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(1,1) = cell2mat(PEPP_MMN_Baseline_DCM(83,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(2,2) = cell2mat(PEPP_MMN_Baseline_DCM(84,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(3,3) = cell2mat(PEPP_MMN_Baseline_DCM(85,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(4,4) = cell2mat(PEPP_MMN_Baseline_DCM(86,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(5,5) = cell2mat(PEPP_MMN_Baseline_DCM(87,avg_fe_pos));
STD_sp2ss_dp2sp_matrix_FE(6,6) = cell2mat(PEPP_MMN_Baseline_DCM(88,avg_fe_pos));

STD_sp2dp_dp2ii_matrix_FE(3,1) = cell2mat(PEPP_MMN_Baseline_DCM(45,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(4,2) = cell2mat(PEPP_MMN_Baseline_DCM(46,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(5,3) = cell2mat(PEPP_MMN_Baseline_DCM(47,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(6,4) = cell2mat(PEPP_MMN_Baseline_DCM(48,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(1,3) = cell2mat(PEPP_MMN_Baseline_DCM(66,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(2,4) = cell2mat(PEPP_MMN_Baseline_DCM(67,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(3,5) = cell2mat(PEPP_MMN_Baseline_DCM(68,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(4,6) = cell2mat(PEPP_MMN_Baseline_DCM(69,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(1,1) = cell2mat(PEPP_MMN_Baseline_DCM(83,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(2,2) = cell2mat(PEPP_MMN_Baseline_DCM(84,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(3,3) = cell2mat(PEPP_MMN_Baseline_DCM(85,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(4,4) = cell2mat(PEPP_MMN_Baseline_DCM(86,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(5,5) = cell2mat(PEPP_MMN_Baseline_DCM(87,avg_fe_pos));
STD_sp2dp_dp2ii_matrix_FE(6,6) = cell2mat(PEPP_MMN_Baseline_DCM(88,avg_fe_pos));

DEV_sp2ss_dp2sp_matrix_FE(3,1) = cell2mat(PEPP_MMN_Baseline_DCM(49,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(4,2) = cell2mat(PEPP_MMN_Baseline_DCM(50,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(5,3) = cell2mat(PEPP_MMN_Baseline_DCM(51,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(6,4) = cell2mat(PEPP_MMN_Baseline_DCM(52,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(1,3) = cell2mat(PEPP_MMN_Baseline_DCM(70,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(2,4) = cell2mat(PEPP_MMN_Baseline_DCM(71,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(3,5) = cell2mat(PEPP_MMN_Baseline_DCM(72,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(4,6) = cell2mat(PEPP_MMN_Baseline_DCM(73,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(1,1) = cell2mat(PEPP_MMN_Baseline_DCM(89,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(2,2) = cell2mat(PEPP_MMN_Baseline_DCM(90,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(3,3) = cell2mat(PEPP_MMN_Baseline_DCM(91,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(4,4) = cell2mat(PEPP_MMN_Baseline_DCM(92,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(5,5) = cell2mat(PEPP_MMN_Baseline_DCM(93,avg_fe_pos));
DEV_sp2ss_dp2sp_matrix_FE(6,6) = cell2mat(PEPP_MMN_Baseline_DCM(94,avg_fe_pos));

DEV_sp2dp_dp2ii_matrix_FE(3,1) = cell2mat(PEPP_MMN_Baseline_DCM(53,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(4,2) = cell2mat(PEPP_MMN_Baseline_DCM(54,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(5,3) = cell2mat(PEPP_MMN_Baseline_DCM(55,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(6,4) = cell2mat(PEPP_MMN_Baseline_DCM(56,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(1,3) = cell2mat(PEPP_MMN_Baseline_DCM(74,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(2,4) = cell2mat(PEPP_MMN_Baseline_DCM(75,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(3,5) = cell2mat(PEPP_MMN_Baseline_DCM(76,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(4,6) = cell2mat(PEPP_MMN_Baseline_DCM(77,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(1,1) = cell2mat(PEPP_MMN_Baseline_DCM(89,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(2,2) = cell2mat(PEPP_MMN_Baseline_DCM(90,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(3,3) = cell2mat(PEPP_MMN_Baseline_DCM(91,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(4,4) = cell2mat(PEPP_MMN_Baseline_DCM(92,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(5,5) = cell2mat(PEPP_MMN_Baseline_DCM(93,avg_fe_pos));
DEV_sp2dp_dp2ii_matrix_FE(6,6) = cell2mat(PEPP_MMN_Baseline_DCM(94,avg_fe_pos));

MMN_matrix_FE(3,1) = cell2mat(PEPP_MMN_Baseline_DCM(57,avg_fe_pos));
MMN_matrix_FE(4,2) = cell2mat(PEPP_MMN_Baseline_DCM(58,avg_fe_pos));
MMN_matrix_FE(5,3) = cell2mat(PEPP_MMN_Baseline_DCM(59,avg_fe_pos));
MMN_matrix_FE(6,4) = cell2mat(PEPP_MMN_Baseline_DCM(60,avg_fe_pos));
MMN_matrix_FE(1,3) = cell2mat(PEPP_MMN_Baseline_DCM(78,avg_fe_pos));
MMN_matrix_FE(2,4) = cell2mat(PEPP_MMN_Baseline_DCM(79,avg_fe_pos));
MMN_matrix_FE(3,5) = cell2mat(PEPP_MMN_Baseline_DCM(80,avg_fe_pos));
MMN_matrix_FE(4,6) = cell2mat(PEPP_MMN_Baseline_DCM(81,avg_fe_pos));
MMN_matrix_FE(1,1) = cell2mat(PEPP_MMN_Baseline_DCM(95,avg_fe_pos));
MMN_matrix_FE(2,2) = cell2mat(PEPP_MMN_Baseline_DCM(96,avg_fe_pos));
MMN_matrix_FE(3,3) = cell2mat(PEPP_MMN_Baseline_DCM(97,avg_fe_pos));
MMN_matrix_FE(4,4) = cell2mat(PEPP_MMN_Baseline_DCM(98,avg_fe_pos));
MMN_matrix_FE(5,5) = cell2mat(PEPP_MMN_Baseline_DCM(99,avg_fe_pos));
MMN_matrix_FE(6,6) = cell2mat(PEPP_MMN_Baseline_DCM(100,avg_fe_pos));
end


% Now plot
if 1 == 1 % silly way to compress it
% STD_sp2ss_dp2sp_matrix_C
figure;imagesc(STD_sp2ss_dp2sp_matrix_C);
xticklabels(ROIs) 
set(gca, 'XAxisLocation', 'top')
yticklabels(ROIs) 
xlabel('From') 
ylabel('To') 
title ([modality_data{mode} ' ' title_plot{tymmn} ' sp2ss/dp2sp HC'])  
% colorbar('eastoutside');
colormap jet
caxis([-.2 0.3])

% DEV_sp2ss_dp2sp_matrix_C
figure;imagesc(DEV_sp2ss_dp2sp_matrix_C);
xticklabels(ROIs) 
set(gca, 'XAxisLocation', 'top')
yticklabels(ROIs) 
xlabel('From') 
ylabel('To') 
title ([modality_data{mode} ' ' title_plot{tymmn} ' sp2ss/dp2sp HC'])  
% colorbar('eastoutside');
colormap jet
caxis([-0.2 0.7])

% STD_sp2dp_dp2ii_matrix_C
figure;imagesc(STD_sp2dp_dp2ii_matrix_C);
xticklabels(ROIs) 
set(gca, 'XAxisLocation', 'top')
yticklabels(ROIs) 
xlabel('From') 
ylabel('To') 
title ([modality_data{mode} ' ' title_plot{tymmn} ' STD sp2dp/dp2ii HC'])  
% colorbar('eastoutside');
colormap jet
caxis([-0.6 0.3])

% DEV_sp2dp_dp2ii_matrix_C
figure;imagesc(DEV_sp2dp_dp2ii_matrix_C);
xticklabels(ROIs) 
set(gca, 'XAxisLocation', 'top')
yticklabels(ROIs) 
xlabel('From') 
ylabel('To') 
title ([modality_data{mode} ' ' title_plot{tymmn} ' sp2dp/dp2ii HC'])  
% colorbar('eastoutside');
colormap jet
caxis([-0.2 0.7])

% MMN_matrix_C
figure;imagesc(MMN_matrix_C);
xticklabels(ROIs) 
set(gca, 'XAxisLocation', 'top')
yticklabels(ROIs) 
xlabel('From') 
ylabel('To') 
title ([modality_data{mode} ' MMN HC'])  
% colorbar('eastoutside');
colormap jet
caxis([-0.2 0.8])

% STD_sp2ss_dp2sp_matrix_FE
figure;imagesc(STD_sp2ss_dp2sp_matrix_FE);
xticklabels(ROIs) 
set(gca, 'XAxisLocation', 'top')
yticklabels(ROIs) 
xlabel('From') 
ylabel('To') 
title ([modality_data{mode} ' ' title_plot{tymmn} ' sp2ss/dp2sp FE'])  
colorbar('eastoutside');
colormap jet
caxis([-0.2 0.3])

% DEV_sp2ss_dp2sp_matrix_FE
figure;imagesc(DEV_sp2ss_dp2sp_matrix_FE);
xticklabels(ROIs) 
set(gca, 'XAxisLocation', 'top')
yticklabels(ROIs) 
xlabel('From') 
ylabel('To') 
title ([modality_data{mode} ' ' title_plot{tymmn} ' sp2ss/dp2sp FE'])  
colorbar('eastoutside');
colormap jet
caxis([-0.2 0.7])

% STD_sp2dp_dp2ii_matrix_FE
figure;imagesc(STD_sp2dp_dp2ii_matrix_FE);
xticklabels(ROIs) 
set(gca, 'XAxisLocation', 'top')
yticklabels(ROIs) 
xlabel('From') 
ylabel('To') 
title ([modality_data{mode} ' ' title_plot{tymmn} ' sp2dp/dp2ii FE'])  
colorbar('eastoutside');
colormap jet
caxis([-0.6 0.3])

% DEV_sp2dp_dp2ii_matrix_FE
figure;imagesc(DEV_sp2dp_dp2ii_matrix_FE);
xticklabels(ROIs) 
set(gca, 'XAxisLocation', 'top')
yticklabels(ROIs) 
xlabel('From') 
ylabel('To') 
title ([modality_data{mode} ' ' title_plot{tymmn} ' sp2dp/dp2ii FE'])  
colorbar('eastoutside');
colormap jet
caxis([-0.2 0.7])

% MMN_matrix_FE
figure;imagesc(MMN_matrix_FE);
xticklabels(ROIs) 
set(gca, 'XAxisLocation', 'top')
yticklabels(ROIs) 
xlabel('From') 
ylabel('To') 
title ([modality_data{mode} ' MMN FE'])  
colorbar('eastoutside');
colormap jet
caxis([-0.2 0.8])
end

end
end
end
