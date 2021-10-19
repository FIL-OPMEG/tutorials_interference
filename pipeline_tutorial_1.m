%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First example tutorial for the manuscript:
%
% 'Interference Suppression Techniques for OPM-based
% MEG: Opportunities and Challenges'. Seymour et al., (2021)
%
% MATLAB scripts were written by
% Dr. Robert Seymour, July 2021 - September 2021
% For enquiries, please contact: rob.seymour@ucl.ac.uk
%
% Tested on MATLAB 2019a, Windows
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Data Paths
% Please download all code dependencies from Zenodo:
% https://doi.org/10.5281/zenodo.5541312
% If you download and extract this archive to the directory one level
% above this script, the paths below should work.
%
% Or alternatively, you can download the latest versions from:
% - Fieldtrip:      https://www.fieldtriptoolbox.org/download/
% - analyse_OPMEG:  https://github.com/neurofractal/analyse_OPMEG
% - NR4M*:          https://github.com/FIL-OPMEG/NR4M
% - optitrack*:     https://github.com/FIL-OPMEG/optitrack
%
% * = private repositories. Email rob.seymour@ucl.ac.uk for access

root = fileparts(which(mfilename));
cd(root);

fieldtripDir    = fullfile(root,'..','tutorial_OPM_scripts','fieldtrip-master');
script_dir      = fullfile(root,'..','tutorial_OPM_scripts','analyse_OPMEG');
mocap_func      = fullfile(root,'..','tutorial_OPM_scripts','optitrack');
NR4M_dir        = fullfile(root,'..','tutorial_OPM_scripts','NR4M');

% Add Fieldtrip to path
disp('Adding the required scripts to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add other scripts to path
addpath(genpath(script_dir));
addpath(genpath(NR4M_dir));
addpath(mocap_func);


%% BIDS data directory. If you extract it to the directory one level above
% where this script lives, this will work:
% Download from https://doi.org/10.5281/zenodo.5539414
data_dir        = fullfile(root,'..','tutorial_OPM_data');


%% Specify save directory
save_dir        = fullfile(data_dir,'results','tutorial1');
mkdir(save_dir);
cd(save_dir);


%% Specify path to head movement data
motive_data  = fullfile(data_dir,'sub-001','ses-002','meg','motion',...
    'sub-001_ses-002_task-aef_run-001_eul_world.csv');


%% Load the OPM data
% Read in the raw data using BIDS
disp('Loading data...');
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'aef';
cfg.bids.sub    = '001';
cfg.bids.ses    = '002';
cfg.bids.run    = '001';
rawData         = ft_opm_create(cfg);


%% Resample to 1000Hz
cfg                 = [];
cfg.resamplefs      = 1000;
[rawData]           = ft_resampledata(cfg, rawData);


%% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0.1 150];
cfg.plot            = 'yes';
cfg.plot_chans      = 'yes';
cfg.plot_ci         = 'no';
cfg.plot_legend     = 'no';
cfg.transparency    = 0.3;
[pow freq]          = ft_opm_psd(cfg,rawData);
ylim([5 1e4])
print('PSD_raw','-dpng','-r400');


%% Load in optitrack data
opti_data = csv2mat_sm(motive_data);

% Convert mm to cm
opti_data = optitrack_to_cm(opti_data);
Fs        = 1000;


%% Sync up opti-track and rawData
[MovementDataOut, OPMdataOut] = syncOptitrackAndOPMdata(opti_data,...
    rawData,'TriggerChannelName','FluxZ-A');


%% Select just the OPM data on the head
cfg             = [];
cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
rawData_MEG     = ft_selectdata(cfg,OPMdataOut);


%% Plot the position (translation) data
pos_data = MovementDataOut.rigidbodies.data(:,4:6);

% Make the first point 0 0 0
pos_data = zero_optitrack_data(pos_data);

% Make time array starting at 0
t = [0:1:length(pos_data)-1]/Fs;

coord = {'X','Y','Z'};
cols = [198 61 61; 142 185 57; 100, 100, 250]./255;
coord_name = {'X: Left-Right','Y: Up-Down','Z: Forward-Back'};

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


%% Regress the optitrack data from the OPM data
ref                 = (MovementDataOut.rigidbodies.data(:,1:6));

% Interpolate any any missing data with previous value
ref = fillmissing(ref,'previous');

% Low-pass filter the optitrack data at 2Hz
[ref]               = ft_preproc_lowpassfilter(...
    ref', 1000, 2, 5);
ref                 = ref';

% Regress using 10s overlapping windows
[rawData_MEG_reg]   = regress_motive_OPMdata(rawData_MEG,ref,10);


%% Compare PSD gain before/after regression step
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 10];
cfg.plot            = 'yes';
cfg.plot_chans      = 'yes';
cfg.plot_ci         = 'no';
cfg.plot_legend     = 'no';
cfg.transparency    = 0.3;
[pow freq]          = ft_opm_psd_compare(cfg,rawData_MEG,rawData_MEG_reg);
print('opti_regression','-dpng','-r400');


%% Homogenous Field Correction
% Please contact t.tierney@ucl.ac.uk for this script
[data_out_mfc, M, chan_inds] = ft_denoise_hfc(rawData_MEG_reg);

% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 80];
cfg.plot            = 'yes';
cfg.plot_chans      = 'yes';
cfg.plot_ci         = 'no';
cfg.plot_legend     = 'no';
cfg.transparency    = 0.3;
[pow freq]          = ft_opm_psd_compare(cfg,rawData_MEG_reg,data_out_mfc);
print('HFC','-dpng','-r400');


%% Data inspection
cfg             = [];
cfg.blocksize   = 30;
cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'vertical';
cfg.colorgroups = 'allblack';
ft_databrowser(cfg,data_out_mfc);


%% Remove DS - channel is bad
cfg                 = [];
cfg.channel         = vertcat(data_out_mfc.label,'-DS-RAD','-DS-TAN');
data_out_mfc        = ft_selectdata(cfg,data_out_mfc);


%% Temporal Filtering

% Spectral Interpolation
cfg                     = [];
cfg.channel             = 'all';
cfg.dftfilter           = 'yes';
cfg.dftfreq             = [50 100 106 120 150];
cfg.dftreplace          = 'neighbour';
cfg.dftbandwidth        = [1 1 1 2 1];
cfg.dftneighbourwidth   = [1 1 1 2 1];
data_out_si             = ft_preprocessing(cfg,data_out_mfc);

% Plot PSD
cfg                     = [];
cfg.channel             = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length        = 10;
cfg.method              = 'tim';
cfg.foi                 = [30 125];
cfg.plot                = 'yes';
cfg.plot_chans          = 'yes';
cfg.plot_ci             = 'no';
cfg.plot_legend         = 'no';
cfg.transparency        = 0.3;
[pow freq]              = ft_opm_psd(cfg,data_out_si);
print('spectral_interp','-dpng','-r400');

% HP-filter @2Hz to attenuate low-freq drifts
cfg                     = [];
cfg.hpfilter            = 'yes';
cfg.hpfreq              = 2;
cfg.filtord             = 5;
cfg.hpinstabilityfix    = 'reduce';
data_out_si_hp          = ft_preprocessing(cfg,data_out_si);

% Plot PSD
cfg                     = [];
cfg.channel             = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length        = 10;
cfg.method              = 'tim';
cfg.foi                 = [0 10];
cfg.plot                = 'yes';
cfg.plot_chans          = 'yes';
cfg.plot_ci             = 'no';
cfg.plot_legend         = 'no';
cfg.transparency        = 0.3;
[pow freq]              = ft_opm_psd_compare(cfg,data_out_mfc,data_out_si_hp);
print('filtering1','-dpng','-r400');
cfg.foi                 = [40 130];
[pow freq]              = ft_opm_psd_compare(cfg,data_out_mfc,data_out_si_hp);
print('filtering2','-dpng','-r400');

% Low Pass Filter @40Hz
cfg                 = [];
cfg.lpfilter        = 'yes';
cfg.lpfreq          = 40;
data_out_si_lp_hp   = ft_preprocessing(cfg,data_out_si_hp);

% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 50];
cfg.plot            = 'yes';
cfg.plot_chans      = 'yes';
cfg.plot_ci         = 'no';
cfg.plot_legend     = 'no';
cfg.transparency    = 0.3;
[pow freq]          = ft_opm_psd(cfg,data_out_si_lp_hp);
print('temporal_filter1','-dpng','-r400');


%% Calculate max field change /s after the various pre-processing steps
[max_FC1,num_secs] = max_FC_calc(OPMdataOut);
[max_FC2,num_secs] = max_FC_calc(rawData_MEG_reg);
[max_FC3,num_secs] = max_FC_calc(data_out_mfc);
[max_FC4,num_secs] = max_FC_calc(data_out_si_lp_hp);

figure;
set(gcf,'Position',[200 200 800 400]);
plot(num_secs,mean(max_FC1,2),'r','LineWidth',2); hold on;
plot(num_secs,mean(max_FC2,2),'g','LineWidth',2);
plot(num_secs,mean(max_FC3,2),'b','LineWidth',2);
plot(num_secs,mean(max_FC4,2),'k','LineWidth',2);
set(gca,'FontSize',20);
ylabel({'Max Field ';'Change (pT/s)'},'FontSize',25);
xlabel('Time (s)','FontSize',25);
set(gca, 'YScale', 'log');
print('max_FC','-dpng','-r400');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the bit of code I used for manual artefact rejection.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cfg             = [];
cfg.blocksize   = 10;
cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'vertical';
cfg.colorgroups = 'allblack';
arft            = ft_databrowser(cfg,data_out_si_hp);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For anyone trying to exactly reproduce the analysis performed in the
% paper, please load arft.mat from the tutorials_interference folder,
% which contains the indices of the time(s) marked as containing artefacts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Turn the highlighted data into nans
arft.artfctdef.reject        = 'nan';
data_out_si_lp_hp_arft       = ft_rejectartifact(arft, data_out_si_lp_hp);


%% Epoch the data
% Here I am using a custom trial function which looks for triggers on a
% specific channel (e.g. NI-TRIG) of the unfiltered raw data loaded earlier
cfg                         = [];
cfg.rawData                 = rawData;
cfg.trialdef.trigchan       = 'NI-TRIG';
cfg.trialdef.downsample     = 1000;
cfg.correct_time            = 0.0;
cfg.trialdef.prestim        = 0.2;        % pre-stimulus interval
cfg.trialdef.poststim       = 0.5;        % post-stimulus interval
cfg.trialfun                = 'OPM_trialfun_usemat';
banana                      = ft_definetrial(cfg);

% Correct for optitrack sync-ing
banana.trl(:,1) = banana.trl(:,1)-round(OPMdataOut.time{1}(1)*1000);
banana.trl(:,2) = banana.trl(:,2)-round(OPMdataOut.time{1}(1)*1000);
trl_index       = banana.trl(:,1);

% Redefines the pre-processed data
data_all        = [];
cfg             = [];
data_all{1}     = ft_redefinetrial(banana,data_out_si_lp_hp_arft);
data_all{2}     = ft_redefinetrial(banana,rawData_MEG);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here I am running the ERF analysis twice - once for the cleaned data and
% once for the raw data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for d = [2 1]
    data = data_all{d};

    % Remove trials with any nan data (i.e. has been marked as artefactual)
    trial2keep   = [];
    trial2reject = [];
    count        = 1;
    count2       = 1;

    for t = 1:length(data.trial)
        result = sum(isnan(data.trial{t}(:)));
        if ~result
            trial2keep(count) = t;
            count=count+1;
        else
            trial2reject(count2) = t;
            count2       = count2+1;
        end
    end

    % The bit that actually removes the bad trials
    try
        cfg               = [];
        cfg.trials        = trial2keep;
        data              = ft_selectdata(cfg,data);
    catch
    end


    % Perform timelockanalysis & baseline-correct
    cfg             = [];
    cfg.channel     = 'all';
    avg_all         = ft_timelockanalysis(cfg,data);

    cfg             = [];
    cfg.baseline    = [-0.1 0];
    [avg_all]       = ft_timelockbaseline(cfg, avg_all);

    % Plot in fT (not included in the paper - just for illustration)
    cfg             = [];
    cfg.parameter   = 'avg';
    cfg.linewidth   = 2;
    cfg.colorgroups = 'allblack';
    ft_databrowser(cfg,avg_all);

    % Perform Students t-test
    epoched_dataset = [];

    for i = 1:length(data.trial)
        epoched_dataset(:,:,i) = data.trial{1,i};
    end

    SE = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
    avg_all.t_value = avg_all.avg./SE;

    % Plot t-value
    cd(save_dir);
    figure;
    set(gcf,'Position',[1 1 800 600]);
    plot(avg_all.time,avg_all.t_value,'k','LineWidth',2);
    set(gca,'FontSize',27);
    ylabel('t-value','FontSize',35);
    xlabel('Time (s)','FontSize',35);
    xlim([-0.1 0.4]);
    ylim([-11 11]);
    print(['ERF' num2str(d)],'-dpng','-r400');


    % Create and Plot 2D Layout (Fieldtrip)
    if d == 2
        cfg             = [];
        cfg.output      = 'lay_123.mat';
        cfg.grad        = data.grad;
        cfg.channel     = data.label;
        %cfg.headshape   = mesh;
        cfg.rotate      = 0;
        cfg.center      = 'yes';
        cfg.projection  = 'polar';
        cfg.channel     =  'all';
        cfg.overlap     = 'keep';
        lay_123         = ft_prepare_layout(cfg);
    end

    % Plot Using ft_topoplotER and change the colormap to RdBu
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    cmap = colormap(flipud(brewermap(64,'RdBu'))); % change the colormap

    cfg             = [];
    cfg.parameter   = 't_value';
    cfg.channel     = ft_channelselection_opm('TAN',rawData);
    cfg.layout      = lay_123;
    cfg.xlim        = [0.08 0.11];
    cfg.zlim        = [-6 6];
    cfg.linewidth   = 2;
    cfg.showlabels  = 'yes';
    cfg.colormap    = cmap;
    cfg.comment     = 'no';
    figure; set(gcf,'Position',[1 1 800 800]);
    ft_topoplotER(cfg,avg_all); hold on;
    c               = colorbar;
    c.Location      = 'southoutside';
    c.FontSize      = 30;
    c.Label.String  = 't-value';
    print(['fieldtmap' num2str(d)],'-dpng','-r400');

end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Source Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute covariance matrix
cfg                  = [];
cfg.covariance       = 'yes';
cfg.vartrllength     = 2;
cfg.covariancewindow = [0 0.5];
avg                  = ft_timelockanalysis(cfg,data);


%% Load MRI & pre-computed headmodel
mri             = ft_read_mri(fullfile(data_dir,'sub-001','ses-002',...
    'anat','001.nii')); mri.coordsys = 'neuromag';

% Load Headmodel
load(fullfile(data_dir,'derivatives','fieldtrip_sourcespace',...
    'sub-001','sub-001_desc-headmodel.mat'));


%%  LF for Virtual electrode analysis
cfg             = [];
cfg.template    = mri;
cfg.nonlinear   = 'yes';
norm            = ft_volumenormalise([],mri);

% Auditory Cortex
pos = [-48 -22 4; 48 -22 4];

% Now we warp the MNI coordinates using the nonlinear warping parameters
posback         = ft_warp_apply(norm.params,pos,'sn2individual');
% xyz positions in individual coordinates
pos_grid        = ft_warp_apply(pinv(norm.initial),posback);

figure;ft_plot_sens(data.grad);
ft_plot_mesh(pos_grid,'vertexcolor','r');


%% Prepare Leadfield
cfg                 = [];
cfg.method          = 'lcmv';
cfg.channel         = data.label;
cfg.grid.pos        = pos_grid;
cfg.grid.unit       = 'mm';
cfg.headmodel       = headmodel;
cfg.grad            = data.grad;
cfg.reducerank      = 2; %(default = 3 for EEG, 2 for MEG)
cfg.normalize       = 'yes' ; %Normalise Leadfield: 'yes' for beamformer
cfg.normalizeparam  = 1;
lf_2                = ft_prepare_leadfield(cfg);

% Concat the leadfields
lf_concat = cat(2,lf_2.leadfield{:});
for k = 1:2
    lf_2.leadfield{k} = lf_concat;
end


%% Source Analysis
cfg                    = [];
cfg.channel            = data.label;
cfg.grad               = data.grad;
cfg.method             = 'lcmv';
cfg.grid               = lf_2;
cfg.headmodel          = headmodel;
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.fixedori      = 'yes';
cfg.lcmv.projectnoise  = 'yes';
cfg.lcmv.lambda        = '0.1%';
sourceall              = ft_sourceanalysis(cfg, avg);

% Find filter from the first .pos
filter123 = cat(1,sourceall.avg.filter{1,:});

VE          = [];
VE.label    = {'A1'};
VE.fsample  = data.fsample;
for subs=1:size(data.trial,2)
    % note that this is the non-filtered "raw" data
    VE.time{subs}       = data.time{subs};
    VE.trial{subs}(:,:) = filter123(:,:)*data.trial{subs}(:,:);
end


%% Perform timelockanalysis
cfg             = [];
avg_VE      	= ft_timelockanalysis([],VE);


%% Perform student t-test
epoched_dataset = [];

for i = 1:length(data.trial)
    epoched_dataset(:,:,i) = VE.trial{1,i};
end

SE              = std(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
avg_VE.t_value  = avg_VE.avg./SE;


%% Plot the VE t-values over time
cfg             = [];
cfg.channel     = avg_VE.label;
cfg.parameter   = 't_value';
cfg.baseline    = [-0.1 0];
cfg.showlegend  = 'yes';
cfg.xlim        = [-0.1 0.4];
cfg.linecolor   = 'k';
cfg.linewidth   = 2;
cfg.ylim        = [-11 11];
figure; ft_singleplotER(cfg,avg_VE)
set(gca,'FontSize',18);
xlabel('Time (s)','FontSize',20);
ylabel('t-value','FontSize',20);
title('');
print('beamforming','-dpng','-r400');
