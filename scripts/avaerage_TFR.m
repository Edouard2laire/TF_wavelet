%% Settings 
addpath('/Users/edelaire1/Documents/software/brainstorm3'); 
if ~brainstorm('status')
    brainstorm
end

close all;

path = '/Users/edelaire1/Documents/Project/CIHR/CIHR_march_2024/TF/PA03';
folder_out = fullfile(path,'figures');

if ~exist(folder_out)
    mkdir(folder_out);
end

bst_plugin('Load', 'TFNIRS', 1);

%% Inputs 

sub     = 'PA03';
mod     = 'nirs';
clus    = 'frontalleft_hem-HbO_res-4';
hem     = 'HbO';

new_frequency = logspace( log10(0.002), log10(0.5), 500);

%% options
options = load(fullfile(path,'options.mat')).options;
options.wavelet = rmfield(options.wavelet,'freqWindow');
options.colormap = 'jet';
options.clim = [0 0.25];
options.wavelet.display.fontscale = 20;

if contains(clus,'frontalleft' ) 
    options.title_tf = sprintf('%s - %s', 'Frontal Left', hem);
elseif contains(clus,'frontalright' ) 
    options.title_tf = sprintf('%s - %s', 'Frontal Right', hem);
else
    options.title_tf = sprintf('%s - %s', clus, hem);
end

sleep_stage  = {'wake';'N1';'N2';'N3';'REM'};
epi_activity = {'bursts', 'spikes_LR', 'spikes_RL', 'spikes_bilat', 'single_s'};

%files = file_find(path, sprintf('%s*mod-%s*clus-%s_hem-%s.mat', sub,  mod,clus, hem), 2, 0 );
files = file_find(path, sprintf('%s*mod-%s*clus-%s.mat', sub,  mod,clus), 2, 0 );

fprintf(' %d files detected \n', length(files));

segments = [];
for iFile = 1:length(files)
    sData = load(files{iFile}).out_data;

    sleep_events    = sData.events( cellfun(@(x) any(strcmp(x,sleep_stage)), {sData.events.label}));
    epi_events      = sData.events( cellfun(@(x) any(strcmp(x,epi_activity)),{sData.events.label}));
    motion_events   = sData.events( contains({sData.events.label}, 'motion'));

    % Extend epi event to account for hemodynamic response 
    for iEvent = 1:length(epi_events)
        epi_events(iEvent)  =  process_ft_wavelet('extendEvent', epi_events(iEvent), 30, 30  );
    end
    
    [sData.WDdata_avg, sData.time] = process_ft_wavelet('removeZero', sData.WDdata_avg,  sData.time );
    process_ft_wavelet('displayTF_Plane',sData.WDdata_avg, sData.time, struct_copy_fields(options,sData.OPTIONS));

    splitting_events = [sleep_events,epi_events, motion_events]; 

    file_segment = process_ft_wavelet('exctractSegment',sData.WDdata_avg,sData.time, splitting_events , sData.events, sData.OPTIONS.wavelet.freqs_analyzed );
    for iSegment = 1:length(file_segment)
        file_segment(iSegment).nAvg(1) = length(sData.cluster.Sensors);
        file_segment(iSegment).nAvg(2:end) = 0;
    end
    segments = [segments ; file_segment ];
end


% View segment i 
iSegment = 15;
options.wavelet.freqs_analyzed = segments(iSegment).freq;
process_ft_wavelet('displayTF_Plane',segments(iSegment).WData,segments(iSegment).time, options);




%% Average within each segment -  different ammount of averaging
% Filter segment less than 60 seconds long
selected_segments = segments([segments.duration] > 90, :);


averaged_segments  = process_ft_wavelet('averageWithinSegment',selected_segments);
resampled_segments = process_ft_wavelet('resampleFrequency',averaged_segments, new_frequency);


averaged_segments = repmat(resampled_segments(1), 1,length(sleep_stage)) ;
disp(' - - - - - - - - - - - -')
fprintf(' %d segments analyed  \n', length(resampled_segments));
for iStage = 1:length(sleep_stage)
    segments_stage = resampled_segments( strcmp({resampled_segments.label}, sleep_stage{iStage}), : );

    fprintf('%s : %d segment. Duration [min, median, max]: %.f s, %.1fs, %.1fs . Total: %.1f s\n',  ...
                sleep_stage{iStage}, ...
                length(segments_stage), ...
                min([segments_stage.duration]), median([segments_stage.duration]), max([segments_stage.duration]),sum([segments_stage.duration]) );

end

disp(' - - - - - - - - - - - -')


for iStage = 1:length(sleep_stage)
    segments_stage = resampled_segments( strcmp({resampled_segments.label}, sleep_stage{iStage}), : );
    averaged_segments(iStage) = process_ft_wavelet('averageBetweenSegment', segments_stage);
end

nAvg = vertcat(averaged_segments.nAvg);
std_err = cat(2,averaged_segments.WDataStd) ./ repmat( sqrt(nAvg(:,3))',  length(new_frequency), 1 );

fig = process_ft_wavelet('displayPowerSpectrum',cat(2,averaged_segments.WData)', ...
                                          std_err', ...
                                          {averaged_segments.label} , ...
                                          averaged_segments(1).freq , ...
                                          options);
saveas(fig,fullfile(folder_out, sprintf('subject-%s_clus-%s_desc-%s.png', sub ,clus, ...
                    'sleepSpectrum')));

%% Average within each segment -  same ammount of averaging

selected_segments   = segments(cellfun(@(x) any(strcmp(sleep_stage, x)), {segments.label}) & ...
                                [segments.duration] > 90, :);

epoched_segments    = process_ft_wavelet('epochSegment',selected_segments, 60, 30);

averaged_segments  = process_ft_wavelet('averageWithinSegment',epoched_segments);
resampled_segments = process_ft_wavelet('resampleFrequency',averaged_segments, new_frequency);


disp(' - - - - - - - - - - - -')
fprintf(' %d segments of 60s analyed  \n', length(resampled_segments));
for iStage = 1:length(sleep_stage)
    segments_stage = resampled_segments( strcmp({resampled_segments.label}, sleep_stage{iStage}), : );
    fprintf('%s : %d segment. \n',  ...
                sleep_stage{iStage}, ...
                length(segments_stage));
end
disp(' - - - - - - - - - - - -')

n_boot = 100;
averaged_segments = repmat(resampled_segments(1), 1,length(sleep_stage)) ;
for iStage = 1:length(sleep_stage)
    segments_stage = resampled_segments( strcmp({resampled_segments.label}, sleep_stage{iStage}), : );
    if length(segments_stage) < 5
        continue;
    end

    sub_sample = repmat(resampled_segments(1), 1,n_boot) ;
    for iboot = 1:n_boot
        idx = randsample(length(segments_stage),5);
        sub_sample(iboot) = process_ft_wavelet('averageBetweenSegment', segments_stage(idx));
    end
    averaged_segments(iStage) = process_ft_wavelet('averageBetweenSegment', sub_sample);
end

fig = process_ft_wavelet('displayPowerSpectrum',cat(2,averaged_segments.WData)', ...
                                          cat(2,averaged_segments.WDataStd)', ...
                                          {averaged_segments.label} , ...
                                          averaged_segments(1).freq , ...
                                          options);

xticks([0.01, 0.1]);
xticklabels({'0.01', '0.1'});

saveas(fig,fullfile(folder_out, sprintf('subject-%s_clus-%s_desc-%s.png', sub ,clus, ...
                    'sleepSpectrum_boot')));

%% burst analysis 

selected_bursts = segments( strcmp({segments.label}, 'bursts') & ...
                                   [segments.duration] == 60, :);

resampled_bursts = process_ft_wavelet('resampleFrequency',selected_bursts,...
                                                         new_frequency);

averaged_bursts  = process_ft_wavelet('averageBetweenSegment', resampled_bursts);

averaged_bursts  = process_ft_wavelet('setTimeOrigin', averaged_bursts, 30 );

viz_option = options;
viz_option.wavelet.freqs_analyzed  = averaged_bursts.freq;
viz_option.title_tf = 'Burst';
viz_option.clim     = [0 0.03]; 

fig = process_ft_wavelet('displayTF_Plane',averaged_bursts.WData,averaged_bursts.time, viz_option);
xline(0, 'r--');
colormap jet; 
clim([0, 0.07])
saveas(fig,fullfile(folder_out, sprintf('subject-%s_clus-%s_desc-%s.png', sub ,clus, ...
                    'TFburst')));

viz_option = options;
viz_option.wavelet.freqs_analyzed  = averaged_bursts.freq;
viz_option.title_tf = 'Burst (z-score)';
viz_option.clim     = [0 3];

zscore = averaged_bursts.WData ./averaged_bursts.WDataStd;
zscore( zscore < 1.96 | zscore > 3 ) = 0;

fig  = process_ft_wavelet('displayTF_Plane',zscore,averaged_bursts.time, viz_option);
colormap jet; 
saveas(fig,fullfile(folder_out, sprintf('subject-%s_clus-%s_desc-%s.png', sub ,clus, ...
                    'TFburst_zscore_thrs')));

save(sprintf('subject-%s_clus-%s_desc-%s.mat', sub ,clus, ...
                    'TFburst'),"averaged_bursts");

