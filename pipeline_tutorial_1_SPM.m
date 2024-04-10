%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First example tutorial for the manuscript:
%
% SPM12 Version
%
% 'Interference Suppression Techniques for OPM-based
% MEG: Opportunities and Challenges'. Seymour et al., (2022). Neuroimage.
%
% MATLAB scripts were written by
% Dr. Robert Seymour, May 2023.
%
% Updated April 2024, using the development version of SPM 
% (downloaded 10/04/2024).
%
% For enquiries, please contact: rob.seymour@ucl.ac.uk
%
% Tested on MATLAB 2019a, Windows
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Code Paths
% This tutorial requires code from:
% - SPM12:          https://www.fil.ion.ucl.ac.uk/spm/software/download/
% - OPM             https://github.com/tierneytim/OPM
% - optitrack:      https://github.com/FIL-OPMEG/optitrack
%
% WORK IN PROGRESS: All code dependencies will soon be packaged up into one
% easy to download .zip file

root = fileparts(matlab.desktop.editor.getActiveFilename);
cd(root);

spm_path        = 'D:\scripts\spm12';
% opm_path        = 'D:\Github\OPM';
mocap_func      = 'D:\Github\optitrack';

% Add SPM12, OPM and optitrack scripts to your path
addpath(spm_path);
addpath(mocap_func);
% addpath(opm_path);

spm_jobman('initcfg');
spm('defaults', 'eeg');

disp('Paths all setup!');


%% BIDS data directory. If you extract it to the directory one level above
% where this script lives, this will work:
% Download from https://doi.org/10.5281/zenodo.5539414
data_dir        = fullfile(root,'..','tutorial_OPM_data_SPM');

%% Specify save directory
save_dir        = fullfile(data_dir,'results','tutorial1');
mkdir(save_dir);
cd(save_dir);

%% Specify path to head movement data
motive_data  = fullfile(data_dir,'sub-001','ses-002','meg','motion',...
    'sub-001_ses-002_task-aef_run-001_eul_world.csv');

%% Specify path to MRI
mri = fullfile(data_dir,'sub-001','ses-002','anat',...
    '001.nii');

%% Load in data
S           = [];
S.data      = fullfile(data_dir,'sub-001','ses-002','meg',...
    'sub-001_ses-002_task-aef_run-001_meg.bin');
S.channels  = fullfile(data_dir,'sub-001','ses-002','meg',...
    'sub-001_ses-002_task-aef_run-001_channels.tsv');
S.meg       = fullfile(data_dir,'sub-001','ses-002','meg',...
    'sub-001_ses-002_task-aef_run-001_meg.json');
S.positions = fullfile(data_dir,'sub-001','ses-002','meg',...
    'sub-001_ses-002_task-aef_run-001_positions.tsv');
S.sMRI      = mri;
D           = spm_opm_create(S)

%% Plot the Data
spm_eeg_review(D)

%% Plot PSD
S             = [];
S.D           = D;
S.triallength = 3000; % time in ms: longer windows provide more frequency resolution but are noisier
S.plot        = 1;
S.channels    = 'MEG';   % select only MEG channels
spm_opm_psd(S);
ylim([1,1e5])         % set y axis limits between 1fT and 100,000 fT
xlim([0.1 140])       % set x axis limits between 0.1 and 140 Hz 
print('PSD_raw','-dpng','-r400');

%% Load in optitrack data
opti_data = csv2mat_sm(motive_data);

% Convert mm to cm
opti_data = optitrack_to_cm(opti_data);
Fs        = 6000;

%% Sync up opti-track and rawData
[MovementDataOut, D_sync] = syncOptitrackAndOPMdata(opti_data,...
    D,'TriggerChannelName',{'FluxZ-A'});

%% Plot the position (translation) data
pos_data = MovementDataOut.rigidbodies.data(:,4:6);

% Make the first point 0 0 0
pos_data = zero_optitrack_data(pos_data);

% Make time array starting at 0
t = [0:1:length(pos_data)-1]/Fs;

coord = {'X','Y','Z'};
cols = [198 61 61; 142 185 57; 100, 100, 250]./255;
coord_name = {'X: Left-Right','Y: Up-Down','Z: Forward-Back'};

% Figure
figure;
set(gcf,'Position',[300 200 1200 500]);

for c = 1:size(pos_data,2)
    m = pos_data(:,c);
    p1 = plot(t,m,'LineWidth',2,'Color',cols(c,:)); hold on;
    p1.Color(4) = 0.7;
end
ylabel(['Distance (cm)'],'FontSize',24);
ylim([-100 100]);
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',24);
legend(coord_name,'Location','EastOutside');
print('opti_pos','-dpng','-r400');

%% Plot the Rotation Data
% Get the rotation data
degXYZ = (MovementDataOut.rigidbodies.data(:,1:3));

% Make the first point 0 0 0
degXYZ = zero_optitrack_data(degXYZ);

% Colours for plotting
cols = [230 0 99;235 210 0; 23 196 230]/255;

coord = {'X','Y','Z'};
coord_name = {'X: Pitch','Y: Roll','Z: Yaw'};

figure;
set(gcf,'Position',[300 200 1200 500]);

for c = 1:size(degXYZ,2)
    m = degXYZ(:,c);
    %m(log_array,:) = NaN;
    p1 = plot(t,m,'LineWidth',2,'Color',cols(c,:)); hold on;
    p1.Color(4) = 0.7;
end
%plot(t,m,'r','LineWidth',2);
ylabel(['Degrees (°)'],'FontSize',24);
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',24);
legend(coord_name,'Location','EastOutside');
print('opti_rot','-dpng','-r400');


%% Remove DS - channel is bad
bad_chans       = {'G2-DS-RAD','G2-DS-TAN'}
list_bad_chans  = find(ismember(D_sync.chanlabels ,bad_chans));
D_sync          = badchannels(D_sync,list_bad_chans,1);

D_sync.save();


%% Homogenous Field Correction
S       = [];
S.D     = D_sync;
S.L     = 1;
[hfD]   = spm_opm_hfc(S);

% Find indices of bad channels
cinds = setdiff(indchantype(hfD,'MEG'),badchannels(hfD));
chans = chanlabels(hfD,cinds);

% Plot relative PSD
S               = [];
S.triallength   = 3000; 
S.plot          = 1;
S.D2            = hfD;
S.D1            = D_sync;
S.channels      = chans;
[shield,freq]   = spm_opm_rpsd(S);
xlim([0.1 140]); hold on;
plot(freq,mean(shield,2),'k','LineWidth',2)
print('HFC','-dpng','-r400');


%% Filter
% HP-filter
S       = [];
S.D     = hfD;
S.freq  = [3];
S.band  = 'high';
fD      = spm_eeg_ffilter(S);

% LP-filter
S       = [];
S.D     = fD;
S.freq  = [40];
S.band  = 'low';
fD      = spm_eeg_ffilter(S);

% Band-stop filter @ 50Hz
S       = [];
S.D     = fD;
S.freq  = [49 51];
S.band  = 'stop';
fD      = spm_eeg_ffilter(S);


%% Plot PSD
S               = [];
S.triallength   = 3000; 
S.plot          = 1;
S.D             = fD;
S.channels      = chans;
[shield,freq]   = spm_opm_psd(S);
ylim([1,1e4]);
xlim([2 40]); hold on;


%% Epoch
S                   = [];
S.D                 = fD;
S.timewin           = [-200 400];
S.triggerChannels   = {'NI-TRIG'};
eD                  = spm_opm_epoch_trigger(S);

eD.save


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the bit of code I used for manual artefact rejection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S       = [];
S.D     = eD;
eD      = spm_eeg_ft_artefact_visual(S);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For anyone trying to exactly reproduce the analysis performed in the
% online tutorial please use these trial indices for badtrials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eD = badtrials(eD,[23 24 44 115 159 185 294 483],1);


%% Average using SPM's robust averaging tool
S           =  [];
S.D         = eD;
muD         = spm_eeg_average(S);

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the ERF - note the peaks at around 100ms corresponding to the 
% M100 auditory evoked field.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Don't plot the bad channels
MEGind  = indchantype(eD,'MEGMAG');
used    = setdiff(MEGind,badchannels(muD));
pl      = muD(used,:,:)';

% Plot the evoked response in fT
figure();
plot(muD.time(),pl,'k','LineWidth',2)
grid on
ax = gca; % current axes
ax.FontSize = 13;
ax.TickLength = [0.02 0.02];
fig= gcf;
fig.Color=[1,1,1];
xlim([-0.1,.30]); % -0.1s to 0.3s
ylim([-220 220]); % -220fT to 220fT
xlabel('Time (s)','FontSize',18)
ylabel('B (fT)','FontSize',18);
print('ERF','-dpng','-r400')

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot a topolpot of the M100 response.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Topoplot (1)
S       = [];
S.D     = muD;
S.T     = .085; 
spm_opm_plotScalpData(S);

%% Topoplot (2)
% Fudge the layout using Fieldtrip
evoked_FT           = muD.ftraw;
evoked_FT.avg       = evoked_FT.trial{1};
evoked_FT.time      = evoked_FT.time{1};
evoked_FT.dimord    = 'chan_time';
evoked_FT           = rmfield(evoked_FT,'trial');
evoked_FT           = rmfield(evoked_FT,'trialinfo');

% Create layout
cfg             = [];
cfg.output      = 'layout.mat';
cfg.grad        = evoked_FT.grad;
cfg.channel     = {'*-TAN'}
cfg.rotate      = 0;
cfg.center      = 'yes';
cfg.projection  = 'polar';
%cfg.overlap     = 'no';
layout        = ft_prepare_layout(cfg);

figure; ft_plot_layout(layout,'fontsize',7)

% Plot Using ft_topoplotER
cfg             = [];
cfg.parameter   = 'avg';
cfg.layout      = lay;
cfg.xlim        = [0.08 0.11];
%cfg.zlim        = [-6 6];
cfg.linewidth   = 2;
cfg.showlabels  = 'yes';
cfg.comment     = 'no';
figure; set(gcf,'Position',[1 1 800 800]);
ft_topoplotER(cfg,evoked_FT); hold on;
c               = colorbar;
c.Location      = 'southoutside';
c.FontSize      = 30;
c.Label.String  = 'fT';
print('topoplot','-dpng','-r400');

%% Topoplot (3) - ANAT

% Get layout from SPM (ANAT Method - Alexander et al., in prep)
fid = fiducials(D);
grad = D.sensors('MEG');
fid_struct = struct('NAS', fid.fid.pnt(contains(fid.fid.label, 'nas'),:), ...
    'LPA', fid.fid.pnt(contains(fid.fid.label, 'lpa'),:), ...
    'RPA', fid.fid.pnt(contains(fid.fid.label, 'rpa'),:));
pos = grad.coilpos;
lay = spm_get_anatomical_layout(pos, grad.label, double(gifti(D.inv{1}.mesh.tess_scalp).vertices), fid_struct, 0);

% Convert to mm
lay = ft_convert_units(lay,'mm');

% Make the mask = to the outline
lay.mask{1} = lay.outline{1};

% Plot Using ft_topoplotER
cfg             = [];
cfg.parameter   = 'avg';
cfg.channel     = {'*-TAN'}
cfg.layout      = lay;
cfg.xlim        = [0.08 0.11];
%cfg.zlim        = [-6 6];
cfg.linewidth   = 2;
cfg.showlabels  = 'yes';
cfg.comment     = 'no';
figure; set(gcf,'Position',[1 1 800 800]);
ft_topoplotER(cfg,evoked_FT); hold on;
c               = colorbar;
c.Location      = 'southoutside';
c.FontSize      = 30;
c.Label.String  = 'fT';
print('topoplot_ANAT','-dpng','-r400');

