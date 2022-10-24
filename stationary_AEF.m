%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analysis pipeline for OPM data collected at WCHN OPM Lab, London
% The data features two participants sitting still, listening to 
% auditory tones.
%
% MATLAB scripts were written by Dr. Robert Seymour, October 2022
% For enquiries, please contact: rob.seymour@ucl.ac.uk
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start-up: Add the appropriate toolboxes and scripts to your path
%           These paths are specific to my PC - change accordingly.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fieldtripDir    = 'D:\scripts\fieldtrip-master';
script_dir      = 'D:\Github\analyse_OPMEG';
mocap_func      = 'D:\Github\optitrack';
atlas_dir       = 'D:\Github\analyse_OPMEG\atlas\HCPMMP';

% Add Fieldtrip to path
disp('Adding Fieldtrip and analyse_OPMEG to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add analyse_OPMEG Scripts to path
addpath(genpath(script_dir));
addpath(mocap_func);

% BIDS data directory. This is specific to my PC - change accordingly.
data_dir        = 'D:\data\auditory_OPM_stationary';
cd(data_dir);

%% Set Subject Number
 
sub = 'sub-002';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_dir        = fullfile(data_dir,'derivatives',sub,'ses-001','results');

motive_data = fullfile(data_dir,sub,'ses-001','mot',...
    [sub '_ses-001_task-aef_run-001_eul_world.csv']);

%% Start preprocessing.
% Read in the raw data using BIDS
disp('Loading data...');
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'aef';
cfg.bids.sub    = sub(end-2:end);
cfg.bids.ses    = '001';
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
cfg.foi             = [0 120];
cfg.plot            = 'yes';
[pow freq]          = ft_opm_psd(cfg,rawData);
ylim([1 1e4])
title('Raw Data');

%% Load in optitrack data
opti_data = csv2mat_sm(motive_data);

% Convert mm to cm
opti_data = optitrack_to_cm(opti_data);

% Plot the rotations
plot_motive_rotation(opti_data,'euler')

% Plot the translations
plot_motive_translation(opti_data,'euler')

% Plot mean marker error (Column 7)
figure;
plot(opti_data.time,opti_data.rigidbodies.data(:,7),'LineWidth',2);
ylabel('Mean Marker Error');xlabel('Time (s)');
drawnow;

%% Sync up opti-track and rawData
[MovementDataOut, OPMdataOut] = syncOptitrackAndOPMdata(opti_data,...
    rawData,'TriggerChannelName','FluxZ-A');

cd(save_dir);

%% Select the OPM Data
cfg             = [];
cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
rawData_MEG     = ft_selectdata(cfg,OPMdataOut);

%% Regress the optitrack data from the MEG data
ref                 = (MovementDataOut.rigidbodies.data(:,1:6));
ref = fillmissing(ref,'previous');
% LP-filter optitrack data
[ref]          = ft_preproc_lowpassfilter(...
    ref', 1000, 2, 5);
ref            = ref';

% Regress
[rawData_MEG_reg]   = regress_motive_OPMdata(rawData_MEG,ref,10);

%% HFC
% Please contact rob.seymour@ucl.ac.uk for this script
[data_out_mfc] = ft_denoise_hfc([],rawData_MEG_reg);
[data_out_mfc_noreg] = ft_denoise_hfc([],rawData_MEG);

%% Plot data
cfg             = [];
cfg.blocksize   = 30;
%cfg.channel     = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.viewmode    = 'butterfly';
cfg.colorgroups = 'allblack';
ft_databrowser(cfg,rawData_MEG);

%% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [0 5];
cfg.plot            = 'yes';
cfg.plot_legend     = 'no';
[pow freq]          = ft_opm_psd_compare(cfg,rawData_MEG,rawData_MEG_reg);
%ylim([1 1e4])

%% Spectral Interpolation
cfg                     = [];
cfg.channel             = 'all';
cfg.dftfilter           = 'yes';
cfg.dftfreq             = [50 83 100 120 150];
cfg.dftreplace          = 'neighbour';
cfg.dftbandwidth        = [2 2 2 3 2];
cfg.dftneighbourwidth   = [1 2 2 2 2];
data_out_si             = ft_preprocessing(cfg,data_out_mfc);

%% HP-filter
cfg                     = [];
cfg.hpfilter            = 'yes';
cfg.hpfreq              = 2;
cfg.filtord             = 5;
cfg.hpinstabilityfix    = 'reduce';
data_out_si_hp          = ft_preprocessing(cfg,data_out_si);

%% Low Pass Filter
cfg                 = [];
cfg.lpfilter        = 'yes';
cfg.lpfreq          = 40;
data_out_si_lp_hp   = ft_preprocessing(cfg,data_out_si_hp);

%% Remove DS (and 17 for 002) - channels are bad
cfg                 = [];
if strcmp(sub,'sub-001')
    cfg.channel         = vertcat(data_out_si_lp_hp.label,'-DS-TAN','-DS-RAD');
elseif strcmp(sub,'sub-002')
    cfg.channel         = vertcat(data_out_si_lp_hp.label,'-17-TAN','-17-RAD','-DS-RAD','-DS-TAN');
end
data_out_si_lp_hp   = ft_selectdata(cfg,data_out_si_lp_hp);
data_out_si_hp      = ft_selectdata(cfg,data_out_si_hp);

%% Trial def
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

% Correct for optitrack
banana.trl(:,1) = banana.trl(:,1)-round(OPMdataOut.time{1}(1)*1000);
banana.trl(:,2) = banana.trl(:,2)-round(OPMdataOut.time{1}(1)*1000);
trl_index       = banana.trl(:,1);

% Redefines the filtered data
cfg     = [];
data    = ft_redefinetrial(banana,data_out_si_lp_hp);

%% Save data
cd(save_dir);
disp('Saving data...');
save('data.mat','data');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Sensor-Level AEF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Perform timelockanalysis
cfg             = [];
cfg.channel    = 'all';
avg_all         = ft_timelockanalysis(cfg,data);

cfg = [];
cfg.baseline = [-0.1 0];
[avg_all] = ft_timelockbaseline(cfg, avg_all);

% Plot in fT
cfg = [];
%cfg.ylim = [-455 455];
cfg.parameter = 'avg';
cfg.linewidth = 2;
cfg.colorgroups   = 'allblack';
ft_databrowser(cfg,avg_all);

% Convert to t-value
epoched_dataset = [];

for i = 1:length(data.trial)
    epoched_dataset(:,:,i) = data.trial{1,i};
end

SE = nanstd(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
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
ylim([-16 16]);

%% Create and Plot 2D Layout (Fieldtrip)
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

%% Select TAN channels
cfg             = [];
cfg.channel    = ft_channelselection_opm('TAN',rawData);
avg_all         = ft_timelockanalysis(cfg,data);

% Select TAN channels from data
cfg             = [];
cfg.channel    = ft_channelselection_opm('TAN',rawData);
data_TAN        = ft_selectdata(cfg,data);

epoched_dataset = [];

for i = 1:length(data_TAN.trial)
    epoched_dataset(:,:,i) = data_TAN.trial{1,i};
end

SE = nanstd(epoched_dataset,[],3)/sqrt(size(epoched_dataset,3));
avg_all.t_value = avg_all.avg./SE;

%% Plot Using ft_topoplotER
% Change the colormap to RdBu
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
cmap = colormap(flipud(brewermap(64,'RdBu'))); % change the colormap

cfg             = [];
cfg.parameter   = 't_value';
cfg.layout      = lay_123;
%cfg.baseline    = [-0.1 0];
cfg.xlim        = [0.08 0.12];
cfg.zlim        = [-8 8];
cfg.linewidth   = 2;
%cfg.zlim        = [-200 200];
cfg.showlabels   = 'yes';
cfg.colormap    = cmap;
cfg.comment     = 'no';
figure; set(gcf,'Position',[1 1 800 800]);
ft_topoplotER(cfg,avg_all); hold on;
c = colorbar;
c.Location = 'southoutside';
c.FontSize = 30;
c.Label.String = 't-value';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Source Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load pre-computed headmodel and sourcemodel
mri             = ft_read_mri(fullfile(data_dir,sub,'ses-001',...
    'anat',[sub(end-2:end) '.nii'])); mri.coordsys = 'neuromag';

scannercast_dir = fullfile(data_dir,'fieldtrip_sourcespace','sub-001')
load(fullfile(data_dir,'derivatives',sub,'ses-001','sourcespace',...
    [sub '_desc-headmodel.mat']));
load(fullfile(data_dir,'derivatives',sub,'ses-001','sourcespace',...
    [sub '_desc-sourcemodel_5mm.mat']));


%% Prepare Leadfield
cfg                 = [];
cfg.method          = 'lcmv';
cfg.channel         = data.label;
cfg.grid            = sourcemodel;
cfg.grid.unit       = 'mm';
cfg.headmodel       = headmodel;
cfg.grad            = data.grad;
cfg.reducerank      = 2%(default = 3 for EEG, 2 for MEG)
cfg.normalize       = 'yes' ; %Normalise Leadfield: 'yes' for beamformer
cfg.normalizeparam  = 1;
lf                  = ft_prepare_leadfield(cfg);

%% Make leadfields symmetric across hemispheres
lf1 = reshape(lf.leadfield, lf.dim);
lf2 = flip(lf1,1);
for k = 1:numel(lf1)
    if ~isempty(lf1{k})&&~isempty(lf2{k})
        lf.leadfield{k} = [lf1{k} lf2{k}];
    else
        lf.leadfield{k} = [];
        lf.inside(k)    = false;
    end
end
clear lf1 lf2

% 
% make a figure of the single subject{i} headmodel, and grid positions
figure; hold on;
ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');
alpha 0.5; camlight;
ft_plot_mesh(lf.pos(lf.inside,:),'vertexsize',1,'vertexcolor','r');
%ft_plot_sens(rawData_MEG.grad, 'style', 'r*'); view([0,0]);
ft_plot_sens(data.grad, 'style', 'r*'); view([0,0]);

% %%
% % the data consists of fewer channels than the precomputed
% % leadfields, the following chunk of code takes care of this
% [a,b] = match_str(data.label, lf.label);
% for k = 1:numel(lf.leadfield)
%     if ~isempty(lf.leadfield{k})
%         tmp = lf.leadfield{k};
%         tmp = tmp(b,:);
%         tmp = tmp-repmat(mean(tmp,1),[size(tmp,1) 1]); % average re-ref
%         lf.leadfield{k} = tmp;
%     end
% end
% lf.label = lf.label(b);

%% Compute covariance matrix
cfg                  = [];
cfg.covariance       = 'yes';
cfg.vartrllength     = 2;
cfg.covariancewindow = [0 0.5];
avg                  = ft_timelockanalysis(cfg,data);

%% Source Analysis
cfg                    = [];
cfg.channel            = data.label;
cfg.grad               = data.grad;
cfg.method             = 'lcmv';
cfg.grid               = lf;
cfg.headmodel          = headmodel;
cfg.lcmv.keepfilter    = 'yes';
cfg.lcmv.fixedori      = 'yes';
cfg.lcmv.projectnoise  = 'yes';
cfg.lcmv.weightnorm    = 'nai';
cfg.lcmv.lambda        = '0.1%';
sourceall              = ft_sourceanalysis(cfg, avg);

% % Remove extra .mom
% for k = 1:size(sourceall.pos,1)
%     if ~isempty(sourceall.avg.mom{k})
%         sourceall.avg.mom{k}(4:6,:) = [];
%     end
% end

%%
% Replace .pos field with template_grid.pos
[t, r] = ft_version;
ddd = load(fullfile(r,'template/sourcemodel/standard_sourcemodel3d5mm.mat'));
template_grid = ddd.sourcemodel;
clear ddd
template_grid = ft_convert_units(template_grid,'mm');

sourceall.pos = template_grid.pos;

%%
% Remove cfg field to save memory
sourceall = rmfield(sourceall,'cfg');

%%
ERF_name = {'M100'};
ERF_toi  = [0.08 0.12]

for ERF = 1:length(ERF_name)
    
    source_pow_post = get_source_pow(data,sourceall,[ERF_toi(ERF,1) ERF_toi(ERF,2)]);
    source_pow_pre  = get_source_pow(data,sourceall,[-0.08 -0.04]);
    % %
    % % source_pow_post.avg.pow = source_pow_post.avg.pow./source_pow_post.avg.noise;
    % % source_pow_pre.avg.pow = source_pow_pre.avg.pow./source_pow_pre.avg.noise;
    % %
    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'pow';
    sourceR = ft_math(cfg,source_pow_post,source_pow_pre);
    
    %% Interpolate
    spm_brain = ft_read_mri('D:\scripts\fieldtrip-master\template\anatomy\single_subj_T1.nii');
    cfg              = [];
    cfg.voxelcoord   = 'no';
    cfg.parameter    = 'pow';
    cfg.interpmethod = 'nearest';
    sourceI  = ft_sourceinterpolate(cfg, sourceR, spm_brain);
    
    %%
    % Change the colormap to RdBu
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    cmap = colormap(flipud(brewermap(64,'RdBu'))); % change the colormap
    
    % Mask bits outside the brain
    %sourceI.anat_mask = spm_brain_seg.brain .* double(sourceI.anatomy);
    
    % Plot
    cfg                 = [];
    cfg.funparameter    = 'pow';
    cfg.funcolormap        = cmap;
    cfg.funcolorlim     = 'maxabs';
    %cfg.maskparameter   = 'anat_mask';
    ft_sourceplot(cfg,sourceI);
    title([ERF_name{ERF}]);
    drawnow;
    
    %% Export to nifti formt and use your favourite MRI software to visualise
    cd(save_dir);
    cfg = [];
    cfg.filetype = 'nifti';
    cfg.filename = [ERF_name{ERF}];
    cfg.parameter = 'pow';
    ft_sourcewrite(cfg,sourceI);
end
