%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Second example tutorial for the manuscript:
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
%
% Or alternatively, you can download the latest versions from:
% - Fieldtrip:     https://www.fieldtriptoolbox.org/download/
% - analyse_OPMEG: https://github.com/neurofractal/analyse_OPMEG
% - NR4M*:          https://github.com/FIL-OPMEG/NR4M
%
% * = private repository. Email rob.seymour@ucl.ac.uk for access

fieldtripDir    = 'D:\data\tutorial_OPM_scripts\fieldtrip-master';
script_dir      = 'D:\data\tutorial_OPM_scripts\analyse_OPMEG';
NR4M_dir        = 'D:\data\tutorial_OPM_scripts\NR4M';

% Add Fieldtrip to path
disp('Adding the required scripts to your MATLAB path');
addpath(fieldtripDir)
ft_defaults;

% Add other scripts to path
addpath(genpath(script_dir));
addpath(genpath(NR4M_dir));


%% BIDS data directory. This is specific to my PC - change accordingly.
% Download from: https://doi.org/10.5281/zenodo.5539414
data_dir        = 'D:\data\tutorial_OPM_data';


%% Specify Save Directory
save_dir        = fullfile(data_dir,'results','tutorial2');
mkdir(save_dir);
cd(save_dir);


%% Read in the raw BIDS-organised data
disp('Loading data...');
cfg             = [];
cfg.folder      = data_dir;
cfg.precision   = 'single';
cfg.bids.task   = 'motor';
cfg.bids.sub    = '001';
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
cfg.foi             = [0.1 150];
cfg.plot            = 'yes';
cfg.plot_chans      = 'yes';
cfg.plot_ci         = 'no';
cfg.plot_legend     = 'no';
cfg.transparency    = 0.3;
[pow freq]          = ft_opm_psd(cfg,rawData);
ylim([5 1e4])
print('raw_data_psd','-dpng','-r300');


%% Filter the data
% Spectral Interpolation
cfg                     = [];
cfg.channel             = 'all';
cfg.dftfilter           = 'yes';
cfg.dftfreq             = [21 83 100];
cfg.dftreplace          = 'neighbour';
cfg.dftbandwidth        = [1 1 1];
cfg.dftneighbourwidth   = [1 1 1];
rawData_si              = ft_preprocessing(cfg,rawData);

% High-pass filter at 2Hz to remove low-frequency artefacts
cfg                     = [];
cfg.hpfilter            = 'yes';
cfg.hpfreq              = 2;
rawData_si_hp           = ft_preprocessing(cfg,rawData_si);

% Low Pass Filter
cfg                     = [];
cfg.lpfilter            = 'yes';
cfg.lpfreq              = 80;
rawData_si_hp_lp        = ft_preprocessing(cfg,rawData_si_hp);


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
[pow freq]          = ft_opm_psd(cfg,rawData_si_hp_lp);
ylim([5 1e4])


%% Synthetic Gradiometry Using 100s overlapping windows
cfg                 = [];
% These are the channels on the head
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData),...
    '-N0-TAN','-N4-TAN','-N0-RAD','-N4-RAD');
% These are the reference channels
cfg.refchannel      = ft_channelselection_opm('MEGREF',rawData);
cfg.filter_ref      = [0 20; 20 80];
cfg.derivative      = 'yes';
cfg.return_all      = 'no';
cfg.winsize         = 100;
data_synth_grad     = ft_opm_synth_gradiometer_window(cfg,rawData_si_hp_lp);

% Plot PSD
cfg                 = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [2 80];
cfg.plot            = 'yes';
cfg.plot_chans      = 'yes';
cfg.plot_ci         = 'no';
cfg.plot_legend     = 'no';
cfg.transparency    = 0.3;
[pow freq]          = ft_opm_psd_compare(cfg,rawData_si_hp_lp,data_synth_grad);
print('synth_grad_psd','-dpng','-r300');


%% HFC
[data_out_mfc, M, chan_inds] = ft_denoise_hfc(data_synth_grad);

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
[pow freq]          = ft_opm_psd_compare(cfg,data_synth_grad,data_out_mfc);
print('HFC','-dpng','-r400');


%% Spectral Interpolation for remaining 50Hz line noise
cfg                     = [];
cfg.channel             = 'all';
cfg.dftfilter           = 'yes';
cfg.dftfreq             = [50];
cfg.dftreplace          = 'neighbour';
cfg.dftbandwidth        = [1];
cfg.dftneighbourwidth   = [1];
data_out_mfc       = ft_preprocessing(cfg,data_out_mfc);


%% Make grad structure only for TAN channels
grad_TAN        = [];
grad_TAN.unit   = 'mm';
count           = 1;

for i = 1:length(data_out_mfc.grad.label)
    if strcmp(data_out_mfc.grad.label{i}(end-2:end),'TAN')
        grad_TAN.chanori(count,:) = data_out_mfc.grad.chanori(i,:);
        grad_TAN.chanpos(count,:) = data_out_mfc.grad.chanpos(i,:);
        grad_TAN.chantype{count}  = data_out_mfc.grad.chantype{i};
        grad_TAN.chanunit{count}  = data_out_mfc.grad.chanunit{i};
        grad_TAN.coilori(count,:) = data_out_mfc.grad.coilori(i,:);
        grad_TAN.coilpos(count,:) = data_out_mfc.grad.coilpos(i,:);
        grad_TAN.label{count}     = data_out_mfc.grad.label{i};
        count                     = count+1;
    end
end

% Create and Plot 2D Layout for TAN channels (Fieldtrip)
cfg             = [];
cfg.output      = 'lay_TAN.mat';
cfg.grad        = grad_TAN;
%cfg.headshape   = mesh;
cfg.rotate      = 0;
cfg.center      = 'yes';
cfg.projection  = 'polar';
cfg.channel     =  'all';
%cfg.overlap     = 'no';
lay_TAN         = ft_prepare_layout(cfg);


%% ICA
% Run ICA
disp('About to run ICA using the fastica method')
cfg                 = [];
cfg.channel         = 'all';
cfg.method          = 'fastica';
cfg.numcomponent    = 50;
cfg.feedback        = 'textbar';
cfg.updatesens      = 'no';
cfg.randomseed      = 454;
comp                = ft_componentanalysis(cfg, data_out_mfc);

% Create Fieldtrip data structure with ICA components instead of MEG data
data_ICA            = data_out_mfc;
data_ICA.label      = comp.label;
data_ICA.trial{1}   = comp.trial{1};


%% Build an interactive IC plot viewer
% Calculate PSD
cfg                 = [];
cfg.channel         = 'all';
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [1 210];
cfg.plot            = 'no';
[pow freq]          = ft_opm_psd(cfg,data_ICA);
po                  = nanmean(pow(:,:,:),3);

% Use Brewermap :colors RdBu
ft_hastoolbox('brewermap',1);
colormap123 = colormap(flipud(brewermap(64,'RdBu')));

% Create Figure
S.f = figure;

for c = 1:30
    set(gcf, 'Position',  [300, 0, 1100, 900]);
    %create two pushbttons
    S.pb = uicontrol('style','push',...
        'units','pix',...
        'position',[450 300 200 40],...
        'fontsize',14,...
        'Tag','flip_button',...
        'string','NEXT',...
        'UserData',struct('flip',0),...
        'callback',@pb_call);
    
    % Find limits of y-axis
    [minDistance, indexOfMin] = min(abs(freq-2));
    [minDistance, indexOfMax] = min(abs(freq-100));
    max_lim = max(po(indexOfMin:indexOfMax,c))*1.1;
    min_lim = min(po(indexOfMin:indexOfMax,c))*1.1;

    % Plot topoplot
    cfg = [];
    cfg.component = c;       % specify the component(s) that should be plotted
    cfg.layout    = lay_TAN; % specify the layout file that should be used for plotting
    cfg.comment   = 'no';
    cfg.marker    = 'labels';
    cfg.colorbar  = 'EastOutside';
    cfg.zlim      = 'maxabs';
    cfg.colormap  = colormap123;
    subplot(8,4,[3:4 7:8 11:12 15:16]);ft_topoplotIC(cfg, comp)
    set(gca,'FontSize',16) % Creates an axes and sets its FontSize to 18
    
    % Plot the raw data from 1-11s
    subplot(8,4,[28:32]);
    plot(comp.time{1},comp.trial{1}(c,1:end),'linewidth',2);
    xlim([1 11]);
    set(gca,'FontSize',16) % Creates an axes and sets its FontSize to 18
    xlabel('Time(s)','FontSize',23);
    ylabel({'Magnetic';'Field (T)'},'FontSize',23);
    
        % Plot PSD
    subplot(8,4,[1:2 5:6 9:10 13:14]);
    semilogy(freq,po(:,c),'-k','LineWidth',2);
    xlim([2 100]);
    set(gca,'FontSize',16) % Creates an axes and sets its FontSize to 18
    xlabel('Frequency (Hz)','FontSize',23)
    labY = ['$$PSD (' 'fT' ' \sqrt[-1]{Hz}$$)'];
    ylabel(labY,'interpreter','latex','FontSize',23)
    ylim([min_lim max_lim]);
    %print(['component_PSD' num2str(c)],'-dpng','-r300');
    
    % If component 6 or 10 save a picture
    if ismember(c,[6 10])
        print(['component_' num2str(c)],'-dpng','-r300');
    end

    % Wait for user input
    uiwait(S.f)
    h = findobj('Tag','moveon');
    
    % Clear the Figure for the next sensor
    clf(S.f);
    
end


%% Post-ICA processing:
% The original data can now be reconstructed, excluding specified components
cfg                 = [];
cfg.component       = [6 10]; %these are the components to be removed
data_clean          = ft_rejectcomponent(cfg, comp,data_out_mfc);

% Plot PSD after ICA
cfg                 = [];
cfg.channel         = 'all';
cfg.trial_length    = 10;
cfg.method          = 'tim';
cfg.foi             = [2 80];
cfg.plot            = 'yes';
cfg.plot_chans      = 'yes';
cfg.plot_ci         = 'no';
cfg.plot_legend     = 'no';
cfg.transparency    = 0.3;
[pow freq]          = ft_opm_psd_compare(cfg,data_out_mfc,data_clean);
print('ICA_gain','-dpng','-r300')


%% Calculate max field change /s following various pre-processign steps

[max_FC1,num_secs] = max_FC_calc(rawData);
[max_FC2,num_secs] = max_FC_calc(rawData_si_hp_lp);
[max_FC3,num_secs] = max_FC_calc(data_synth_grad);
[max_FC4,num_secs] = max_FC_calc(data_out_mfc);
[max_FC5,num_secs] = max_FC_calc(data_clean);

% Make Figure to show the results
figure; 
set(gcf,'Position',[200 200 800 400]);
plot(num_secs,mean(max_FC1,2),'r','LineWidth',2); hold on;
plot(num_secs,mean(max_FC2,2),'g','LineWidth',2);
plot(num_secs,mean(max_FC3,2),'b','LineWidth',2);
plot(num_secs,mean(max_FC4,2),'k','LineWidth',2);
plot(num_secs,mean(max_FC5,2),'Color',[0.4940, 0.1840, 0.5560],'LineWidth',2);
set(gca,'FontSize',20);
ylabel({'Max Field ';'Change (pT/s)'},'FontSize',25);
xlabel('Time (s)','FontSize',25);
xlim([2 800]);
set(gca, 'YScale', 'log');
print('maxFC','-dpng','-r300');


%% Epoch the data
cfg                         = [];
cfg.rawData                 = rawData;
cfg.trialdef.trigchan       = 'NI-TRIG';
cfg.trialdef.downsample     = 1000;
cfg.correct_time            = [];
cfg.trialdef.prestim        = 2.0;        % pre-stimulus interval
cfg.trialdef.poststim       = 6.0;        % post-stimulus interval
cfg.trialfun                = 'OPM_trialfun_usemat';
banana                      = ft_definetrial(cfg);

% Redefines the filtered data
cfg         = [];
data        = ft_redefinetrial(banana,data_clean);


%% Calculate TFRs using a hanning taper
cfg              = [];
cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.pad          = 'nextpow2';
cfg.foi          = 1:2:41;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -2:0.05:6;
TFR              = ft_freqanalysis(cfg, data);


%% Plot Fieldmaps for beta desync and beta rebound
cfg                 = [];
cfg.colormap        = colormap123;
cfg.channel         = vertcat(ft_channelselection_opm('TAN',rawData));
cfg.layout          = lay_TAN;
cfg.ylim            = [13 30];
cfg.xlim            = [0 2.5];
cfg.baseline        = [-1.5 0];
cfg.baselinetype    = 'db';
cfg.zlim            = 'maxabs';
cfg.comment         = 'no';
figure;ft_topoplotTFR(cfg,TFR); hold on;
c                   = colorbar;
c.FontSize          = 18;
c.Label.String = 'Power (dB)';
print('beta_desync','-dpng','-r300');
cfg.xlim            = [2.5 4];
figure;ft_topoplotTFR(cfg,TFR); hold on;
c                   = colorbar;
c.FontSize          = 18;
c.Label.String      = 'Power (dB)';
print('beta_rebound','-dpng','-r300');


%% Plot sensors with max dB values for beta desync & rebound
% Sensor: MZ-TAN
cfg                 = [];
cfg.channel         = 'MZ-TAN';
cfg.colormap        = colormap123;
cfg.layout          = lay_TAN;
cfg.ylim            = [1 41];
cfg.xlim            = [-0.5 5.5];
cfg.baseline        = [-1.5 0];
cfg.baselinetype    = 'db';
cfg.zlim            = 'maxabs';
cfg.comment         = 'no';
figure;ft_singleplotTFR(cfg,TFR); hold on;
title('');
c                   = colorbar;
c.FontSize          = 18;
c.Label.String      = 'Power (dB)';
set(gca,'FontSize',20);
xlabel('Time (s)','FontSize',25);
ylabel('Frequency (Hz)','FontSize',25);
plot(repmat(2.5,length([1:2:41])),[1:2:41],'--k','LineWidth',2);
print('beta_desync_MZ','-dpng','-r300');

% Now DQ
cfg.channel         = 'DQ-TAN';
figure;ft_singleplotTFR(cfg,TFR); hold on;
title('');
c                   = colorbar;
c.FontSize          = 18;
c.Label.String      = 'Power (dB)';
set(gca,'FontSize',20);
xlabel('Time (s)','FontSize',25);
ylabel('Frequency (Hz)','FontSize',25);
plot(repmat(2.5,length([1:2:41])),[1:2:41],'--k','LineWidth',2);
print('beta_desync_DQ','-dpng','-r300');


%% Repeat the TFR analysis with raw, unprocessed data
% Redefines the filtered data
cfg              = [];
data_raw         = ft_redefinetrial(banana,rawData);

% TFR
cfg              = [];
cfg.channel         = vertcat(ft_channelselection_opm('MEG',rawData));

cfg.output       = 'pow';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.pad          = 'nextpow2';
cfg.foi          = 1:2:41;
cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5;   % length of time window = 0.5 sec
cfg.toi          = -2:0.05:6;
TFR              = ft_freqanalysis(cfg, data_raw);

% Sensor: MZ-TAN
cfg                 = [];
cfg.channel         = 'MZ-TAN';
cfg.colormap        = colormap123;
cfg.layout          = lay_TAN;
cfg.ylim            = [1 41];
cfg.xlim            = [-0.5 5.5];
cfg.baseline        = [-1.5 0];
cfg.baselinetype    = 'db';
cfg.zlim = [-8 8];
cfg.comment         = 'no';
figure;ft_singleplotTFR(cfg,TFR); hold on;
title('');
c                   = colorbar;
c.FontSize          = 18;
c.Label.String      = 'Power (dB)';
set(gca,'FontSize',20);
xlabel('Time (s)','FontSize',25);
ylabel('Frequency (Hz)','FontSize',25);
plot(repmat(2.5,length([1:2:41])),[1:2:41],'--k','LineWidth',2);
print('beta_desync_MZ_raw','-dpng','-r300');

% Now DQ
cfg.channel         = 'DQ-TAN';
figure;ft_singleplotTFR(cfg,TFR); hold on;
title('');
c                   = colorbar;
c.FontSize          = 18;
c.Label.String      = 'Power (dB)';
set(gca,'FontSize',20);
xlabel('Time (s)','FontSize',25);
ylabel('Frequency (Hz)','FontSize',25);
plot(repmat(2.5,length([1:2:41])),[1:2:41],'--k','LineWidth',2);
print('beta_desync_DQ_raw','-dpng','-r300');

