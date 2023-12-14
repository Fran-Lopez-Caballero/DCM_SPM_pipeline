%% First, manually openfig in SPM each GAVR from each group and save each condition in the path

%% Plot EEG
clear;close all;
% Load subject array and calculate n
load('Private/Path\DCM\Scripts_final_pipeline\subject_array_EEG_matched.mat');
count_C_EEG = 0;
for i = 1:size(subject_array,1)
    if strcmp(subject_array{i,2},'C') && ~strcmp(subject_array{i,5},'bad_EEG')
        count_C_EEG = count_C_EEG + 1;
    end
end
count_FE_EEG = 0;
for i = 1:size(subject_array,1)
    if strcmp(subject_array{i,2},'FE') && ~strcmp(subject_array{i,5},'bad_EEG')
        count_FE_EEG = count_FE_EEG + 1;
    end
end

% Retrieve data and x axis
openfig('Private/Path/DCM/SPM_MEG_data/QC/GAVR_butterworth_filter/C_pMMN_EEG.fig');
l=findobj(gca,'type','line');
time_axis = get(l(1),'xdata');
C_pMMN_EEG =get(l(1),'ydata');
max_so_far = max(C_pMMN_EEG);
min_so_far = min(C_pMMN_EEG);
close all;
openfig('Private/Path\DCM\SPM_MEG_data\QC\GAVR_butterworth_filter\FE_pMMN_EEG.fig');
l=findobj(gca,'type','line');
FE_pMMN_EEG =get(l(1),'ydata');
if max(FE_pMMN_EEG) > max_so_far
    max_so_far = max(FE_pMMN_EEG);
end
if min(FE_pMMN_EEG) < min_so_far
    min_so_far = min(FE_pMMN_EEG);
end
close all;
openfig('Private/Path\DCM\SPM_MEG_data\QC\GAVR_butterworth_filter\C_dMMN_EEG.fig');
l=findobj(gca,'type','line');
C_dMMN_EEG =get(l(1),'ydata');
if max(C_dMMN_EEG) > max_so_far
    max_so_far = max(C_dMMN_EEG);
end
if min(C_dMMN_EEG) < min_so_far
    min_so_far = min(C_dMMN_EEG);
end
close all;
openfig('Private/Path\DCM\SPM_MEG_data\QC\GAVR_butterworth_filter\FE_dMMN_EEG.fig');
l=findobj(gca,'type','line');
FE_dMMN_EEG =get(l(1),'ydata');
if max(FE_dMMN_EEG) > max_so_far
    max_so_far = max(FE_dMMN_EEG);
end
if min(FE_dMMN_EEG) < min_so_far
    min_so_far = min(FE_dMMN_EEG);
end

ylimitmax = max_so_far+(max_so_far*0.05);
ylimitmin = min_so_far+(min_so_far*0.05);

% Plot in a 1x2 figure
set(0,'defaultfigurecolor',[1 1 1]);
close all;
figure;
num_subplots = 2; % For later
h(1) = subplot(1,2,1); h(2) = subplot(1,2,2);
hold (h(1),'on')
plot(h(1),time_axis,C_pMMN_EEG);
hold (h(1),'on')
plot(h(1),time_axis,FE_pMMN_EEG);
hold (h(2),'on')
plot(h(2),time_axis,C_dMMN_EEG);
hold (h(2),'on')
plot(h(2),time_axis,FE_dMMN_EEG);
legend_text = {};
legend_text{1} = ['C (n = ' num2str(count_C_EEG) ')'];
legend_text{2} = ['FE (n = ' num2str(count_FE_EEG) ')'];
for i = 1:num_subplots
    ylim(h(i),[ylimitmin,ylimitmax]);
    xlim(h(i),[-50 300]);
    line(h(i),[-50 300],[0 0], 'Color','black')
    line(h(i),[0 0], [ylimitmin ylimitmax],'Color','black')
end
ylabel(h(1),'field intensity (in µV)'); 
xlabel(h(1),'Time (in ms)'); 
legend(h(1),legend_text);
suptitle('EEG (FCz) SPM analysis')

%% Plot MEG

clear;close all;
% Load subject array and calculate n
load('Private/Path\DCM\Scripts_final_pipeline\subject_array_EEG_matched.mat');
count_C_MEG = 0;
for i = 1:size(subject_array,1)
    if strcmp(subject_array{i,2},'C') && ~strcmp(subject_array{i,6},'bad_MEG')
        count_C_MEG = count_C_MEG + 1;
    end
end
count_FE_MEG = 0;
for i = 1:size(subject_array,1)
    if strcmp(subject_array{i,2},'FE') && ~strcmp(subject_array{i,6},'bad_MEG')
        count_FE_MEG = count_FE_MEG + 1;
    end
end

% Retrieve data and x axis
openfig('Private/Path/DCM/SPM_MEG_data/QC/GAVR_butterworth_filter/C_pMMN_MEG.fig');
l=findobj(gca,'type','line');
time_axis = get(l(1),'xdata');
C_pMMN_MEG =get(l(1),'ydata');
max_so_far = max(C_pMMN_MEG);
min_so_far = min(C_pMMN_MEG);
close all;
openfig('Private/Path\DCM\SPM_MEG_data\QC\GAVR_butterworth_filter\FE_pMMN_MEG.fig');
l=findobj(gca,'type','line');
FE_pMMN_MEG =get(l(1),'ydata');
if max(FE_pMMN_MEG) > max_so_far
    max_so_far = max(FE_pMMN_MEG);
end
if min(FE_pMMN_MEG) < min_so_far
    min_so_far = min(FE_pMMN_MEG);
end
close all;
openfig('Private/Path\DCM\SPM_MEG_data\QC\GAVR_butterworth_filter\C_dMMN_MEG.fig');
l=findobj(gca,'type','line');
C_dMMN_MEG =get(l(1),'ydata');
if max(C_dMMN_MEG) > max_so_far
    max_so_far = max(C_dMMN_MEG);
end
if min(C_dMMN_MEG) < min_so_far
    min_so_far = min(C_dMMN_MEG);
end
close all;
openfig('Private/Path\DCM\SPM_MEG_data\QC\GAVR_butterworth_filter\FE_dMMN_MEG.fig');
l=findobj(gca,'type','line');
FE_dMMN_MEG =get(l(1),'ydata');
if max(FE_dMMN_MEG) > max_so_far
    max_so_far = max(FE_dMMN_MEG);
end
if min(FE_dMMN_MEG) < min_so_far
    min_so_far = min(FE_dMMN_MEG);
end

ylimitmax = max_so_far+(max_so_far*0.05);
ylimitmin = min_so_far+(min_so_far*0.05);

% Plot in a 1x2 figure
set(0,'defaultfigurecolor',[1 1 1]);
close all;
figure;
num_subplots = 2; % For later
h(1) = subplot(1,2,1); h(2) = subplot(1,2,2);
hold (h(1),'on')
plot(h(1),time_axis,C_pMMN_MEG);
hold (h(1),'on')
plot(h(1),time_axis,FE_pMMN_MEG);
hold (h(2),'on')
plot(h(2),time_axis,C_dMMN_MEG);
hold (h(2),'on')
plot(h(2),time_axis,FE_dMMN_MEG);
legend_text = {};
legend_text{1} = ['C (n = ' num2str(count_C_MEG) ')'];
legend_text{2} = ['FE (n = ' num2str(count_FE_MEG) ')'];
for i = 1:num_subplots
    ylim(h(i),[ylimitmin,ylimitmax]);
    xlim(h(i),[-50 300]);
    line(h(i),[-50 300],[0 0], 'Color','black')
    line(h(i),[0 0], [ylimitmin ylimitmax],'Color','black')
end
ylabel(h(1),'field intensity (in fT)'); 
xlabel(h(1),'Time (in ms)'); 
legend(h(1),legend_text);
suptitle('MEG (MEG 2411) SPM analysis')

%% Plot EEG matched DCM

clear;close all;

% Retrieve data and x axis
% STD
openfig('Private/Path/DCM/SPM_MEG_data/Statistics_scalp/EEG_GAVRs_matched/GAVR_EEG_STD_C.fig');
l=findobj(gca,'type','line');
time_axis = get(l(1),'xdata');
C_pMMN_EEG =get(l(1),'ydata');
max_so_far = max(C_pMMN_EEG);
min_so_far = min(C_pMMN_EEG);
close all;
openfig('Private/Path/DCM/SPM_MEG_data/Statistics_scalp/EEG_GAVRs_matched/GAVR_EEG_STD_FE.fig');
l=findobj(gca,'type','line');
FE_pMMN_EEG =get(l(1),'ydata');
if max(FE_pMMN_EEG) > max_so_far
    max_so_far = max(FE_pMMN_EEG);
end
if min(FE_pMMN_EEG) < min_so_far
    min_so_far = min(FE_pMMN_EEG);
end
close all;
openfig('Private/Path/DCM/SPM_MEG_data/Statistics_scalp/EEG_GAVRs_matched/GAVR_EEG_PDev_C.fig');
l=findobj(gca,'type','line');
C_dMMN_EEG =get(l(1),'ydata');
if max(C_dMMN_EEG) > max_so_far
    max_so_far = max(C_dMMN_EEG);
end
if min(C_dMMN_EEG) < min_so_far
    min_so_far = min(C_dMMN_EEG);
end
close all;
openfig('Private/Path/DCM/SPM_MEG_data/Statistics_scalp/EEG_GAVRs_matched/GAVR_EEG_PDev_FE.fig');
l=findobj(gca,'type','line');
FE_dMMN_EEG =get(l(1),'ydata');
if max(FE_dMMN_EEG) > max_so_far
    max_so_far = max(FE_dMMN_EEG);
end
if min(FE_dMMN_EEG) < min_so_far
    min_so_far = min(FE_dMMN_EEG);
end

ylimitmax = .95; % max_so_far+(max_so_far*0.05);
ylimitmin = -1.62; % min_so_far+(min_so_far*0.05);

% Plot in a 1x2 figure
set(0,'defaultfigurecolor',[1 1 1]);
close all;
figure;
num_subplots = 2; % For later
h(1) = subplot(1,2,1); h(2) = subplot(1,2,2);
hold (h(1),'on')
plot(h(1),time_axis,C_pMMN_EEG,'Color',[0.3,0.3,0.3]);
hold (h(1),'on')
plot(h(1),time_axis,FE_pMMN_EEG, 'Color',[0.3,0.3,0.3],'LineStyle','--');
hold (h(2),'on')
plot(h(2),time_axis,C_dMMN_EEG, 'Color',[0.3,0.3,0.3]);
hold (h(2),'on')
plot(h(2),time_axis,FE_dMMN_EEG, 'Color',[0.3,0.3,0.3],'LineStyle','--');
legend_text = {};
legend_text{1} = ['C (n = ' num2str(26) ')'];
legend_text{2} = ['FE (n = ' num2str(26) ')'];
for i = 1:num_subplots
    ylim(h(i),[ylimitmin,ylimitmax]);
    xlim(h(i),[-50 300]);
    line(h(i),[-50 300],[0 0], 'Color','black')
    line(h(i),[0 0], [ylimitmin ylimitmax],'Color','black')
end
ylabel(h(1),'field intensity (in µV)'); 
xlabel(h(1),'Time (in ms)'); 
legend(h(1),legend_text);
suptitle('EEG (FCz) SPM analysis')


gray = [0 0 0];
patch(h(i),[142 142 142 142],[min(ylim)*[1 1] max(ylim)*[1 1]],gray,'FaceAlpha', 1,'EdgeAlpha',0.5,'HandleVisibility','off')
% [-40 0 0 -40]



