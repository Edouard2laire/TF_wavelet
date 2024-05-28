function varargout = process_ft_wavelet( varargin )
% process_ft_wavelet: Call Jean-Marc Continous wavelet transformation for
% time-frequency analysis
%
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% 
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Jean-Marc Lina, Edouard Delaire (2023-2024)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Time-frequency: Morse Wavelet';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = {'NIRS','Frequency'};
    sProcess.Index       = 1705;
    sProcess.Description = 'http://www.fieldtriptoolbox.org/tutorial/timefrequencyanalysis';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import','raw', 'data'};
    sProcess.OutputTypes = {'import','data', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'data', 'raw', 'matrix'};
    sProcess.OutputTypes = {'timefreq', 'timefreq', 'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Options: Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    sProcess.options.timewindow.Group   = 'input';
    % Options: Sensor types

    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'NIRS';

    % === CLUSTERS
    sProcess.options.clusters.Comment = '';
    sProcess.options.clusters.Type    = 'cluster';
    sProcess.options.clusters.Value   = [];


    sProcess.options.label1.Comment = '<b>Wavelet parameters:</b>';
    sProcess.options.label1.Type    = 'label';


    sProcess.options.vanish_moments.Comment = 'Vanishing Moment: ';
    sProcess.options.vanish_moments.Type    = 'value';
    sProcess.options.vanish_moments.Value   = {5.3, '',1};

    sProcess.options.order.Comment = 'Order: ';
    sProcess.options.order.Type    = 'value';
    sProcess.options.order.Value   = {40, '',0};

    sProcess.options.nb_levels.Comment = 'Nb levels: ';
    sProcess.options.nb_levels.Type    = 'value';
    sProcess.options.nb_levels.Value   = {1024, '',0};

    sProcess.options.label2.Comment = '<b>Display options:</b>';
    sProcess.options.label2.Type    = 'label';


    sProcess.options.normalization.Comment = 'Apply normalization: ';
    sProcess.options.normalization.Type    = 'checkbox';
    sProcess.options.normalization.Value   = 1;

    sProcess.options.freq_range.Comment = 'Frequency Range: ';
    sProcess.options.freq_range.Type    = 'range';
    sProcess.options.freq_range.Value   = {[0.002, 0.5], 'Hz', 3};

end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

    % Load recordings
    if strcmp(sInputs.FileType, 'data')     % Imported data structure
        sData = in_bst_data(sInputs(1).FileName);
    elseif strcmp(sInputs.FileType, 'raw')  % Continuous data file       
        sData = in_bst(sInputs(1).FileName, [], 1, 1, 'no');
    end

    % Prepare Data and Options 
    sChannels        = in_bst_channel(sInputs.ChannelFile);
    channelTypes     = sProcess.options.sensortypes.Value;
    freq_range       = sProcess.options.freq_range.Value{1};

    iChannels         = good_channel(sChannels.Channel, sData.ChannelFlag, channelTypes) ;
    
    icluster          = cellfun(@(x)find(strcmp( {sChannels.Clusters.Label},x)),  sProcess.options.clusters.Value);
    cluster           = sChannels.Clusters(icluster);

    if isfield(sProcess.options, 'timewindow') && isfield(sProcess.options.timewindow, 'Value') && iscell(sProcess.options.timewindow.Value) && ~isempty(sProcess.options.timewindow.Value)
        TimeWindow = sProcess.options.timewindow.Value{1};
    end

    if isempty(TimeWindow)
        time = sData.Time;
        iTime = 1:length(time);
    else
        iTime = panel_time('GetTimeIndices', sData.Time, TimeWindow');
        time  = sData.Time(iTime);
    end

    F           = sData.F(iChannels,iTime);
    if strcmpi(channelTypes,'nirs')
        F = nst_misc_convert_to_mumol(F,sData.DisplayUnits);
    else
        F = bst_getunits(F, channelTypes, [], sData.DisplayUnits);
    end

    sChannels   = sChannels.Channel(iChannels);

    p = nextpow2(length(time));
    F = padarray(F,[0 round((2^p - size(F,2))/2)],0,'both');
    ext_time    = padarray(time,[0 round((2^p - size(time,2))/2)],NaN,'both');
    iOrigTime   = ~isnan(ext_time);

    OPTIONS.wavelet.vanish_moments =  sProcess.options.vanish_moments.Value{1} ; % vanish momemts
    OPTIONS.wavelet.order          =  sProcess.options.order.Value{1} ;  % spectral decay
    OPTIONS.wavelet.nb_levels      =  sProcess.options.nb_levels.Value{1};% number of voices
    OPTIONS.wavelet.verbose        = 1;   % verbose or not
    OPTIONS.mandatory.DataTime     = time; 
    OPTIONS.wavelet.display.fontscale = 16;

    if sProcess.options.normalization.Value
        OPTIONS.wavelet.display.TaegerK = 'yes';
    else
        OPTIONS.wavelet.display.TaegerK = 'no';
    end



    % Start of the main code
    wData           = zeros(length(cluster),  OPTIONS.wavelet.nb_levels  + 1,length(time)) ; % N_channel x Nfreq x Ntime
    OPTIONS         = repmat(OPTIONS, 1,length(cluster) );

    nSensors = sum(cellfun(@(x)length(x),{cluster.Sensors}));
    bst_progress('start', 'Running Time-Frequency Analysis', 'Running Time-Frequency Analysis', 0, nSensors);

    for iCluster = 1:length(cluster) 
        wData_temp  = nan(length(cluster(iCluster).Sensors),  OPTIONS(iCluster).wavelet.nb_levels  + 1,length(time)) ; % N_channel x Nfreq x Ntime

        for iSensor = 1:length(cluster(iCluster).Sensors)
            
            % Step 1 - compute time-frequency representation
            iChannel = find(strcmp({sChannels.Name}, cluster(iCluster).Sensors(iSensor)));
            [tmp, OPTIONS(iCluster)] = be_CWavelet(squeeze(F(iChannel,1:2^p)), OPTIONS(iCluster));
            
            % Step 2- Normalize the TF maps ( Remove 1/f)
            [power, title_tf] = normalize(tmp(iOrigTime,:)', OPTIONS(iCluster));
            OPTIONS(iCluster).title_tf =  title_tf ;

            % Step 3- Normalize the TF maps ( standardize power)
            power_time              = sqrt(sum(power.^2));
            wData_temp(iSensor,:,:) = power ./ median(power_time) ;

            bst_progress('inc', 1); 
        end

        % Step 4 - Average accross 
        if length(cluster(iCluster).Sensors) > 1
            wData(iCluster,:,:) = squeeze(mean(wData_temp));
        else
            wData(iCluster,:,:) = squeeze(wData_temp);
        end

        OPTIONS(iCluster).title_tf  = sprintf('%s - %s', strrep(cluster(iCluster).Label,'_',' '), OPTIONS(iCluster).title_tf);
        
        % Step 4. Select frequency Band
        OPTIONS(iCluster) = select_frequency_band(freq_range(1),freq_range(2),OPTIONS(iCluster));

        % Step 5. Visualize the time-frequnecy map
        displayTF_Plane(squeeze(wData(iCluster,:,:)),time, OPTIONS(iCluster));


        % Step 6. Save in Brainstorm
        nfreq = length(OPTIONS(iCluster).wavelet.freqs_analyzed);
        nTime = length(time);
        
        TFmask  = ones(nfreq,nTime);
        is_zero = all(squeeze(wData(iCluster,:,:)) == 0);
        TFmask(:, is_zero) = 0;
         
        sDataOut            = db_template('timefreqmat');
        sDataOut.TF         = permute(wData_temp, [1,3,2]); %nSensor x nTime x nFreq
        sDataOut.RowNames   = cluster(iCluster).Sensors;
        sDataOut.TFmask     = TFmask;
        sDataOut.DataType   = 'data';
        sDataOut.DataFile   = sInputs.FileName;
        sDataOut.Time       = time;
        sDataOut.Freqs      = OPTIONS(iCluster).wavelet.freqs_analyzed;
        sDataOut.Comment    = sprintf('CW [%s]',cluster(iCluster).Label);
        sDataOut.Method     = 'morlet';
        sDataOut.Measure    = 'power';
        sDataOut.Options.PowerUnits = 'physical';
        sDataOut.DisplayUnits       =' \mumol.l-1'; % TODO: fix display units

        % Generate a new file name in the same folder
        sStudy = bst_get('Study', sInputs.iStudy);
        OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'timefreq_wavelet');
        sDataOut.FileName = file_short(OutputFile);
        bst_save(OutputFile, sDataOut, 'v7');
        % Register in database
        db_add_data(sInputs.iStudy, OutputFile, sDataOut);
        OutputFiles{end+1} = OutputFile;

    end 

    bst_progress('stop');
    
end 

function [power, time] = removeZero(power, time )
% Remove the zeros comming from the edge effect 

    is_zero = all(power == 0);
    power = power(:, ~is_zero);
    time  = time(~is_zero);

end

function idx  =  getIndexEvent(Time, Event )
    idx = zeros(1 , length(Time)); 
         
    for i_event = 1:size(Event.times,2)
        i_intra_event = panel_time('GetTimeIndices', Time, Event.times(:,i_event)');
        idx(1,i_intra_event) = 1;
    end

    idx = logical(idx);
end


function Event  =  extendEvent(Event, before, after  )

    if size(Event.times,1) == 1
        Event.times = repmat(Event.times,2,1);
    end
    
    Event.times(1,:) = Event.times(1,:) - before;
    Event.times(2,:) = Event.times(2,:) + after;
    Event.times = Event.times(:, Event.times(2,:) > Event.times(1,:));

end

function [power, title_tf] = normalize(tf, OPTIONS)

    ofs = ceil(0.02*size(tf,2));
    power = zeros(size(tf));

    % On calcule a puissance de Taeger-Kaiser (si demandee):
    if isfield(OPTIONS.wavelet.display,'TaegerK') && strcmp(OPTIONS.wavelet.display.TaegerK,'yes') 
        tf1 = 0.25*tf(:,3:end).*conj(tf(:,1:end-2))+0.25*conj(tf(:,3:end)).*tf(:,1:end-2);       
        power(:,2+ofs:end-1-ofs) = 0.5*abs(tf(:,2+ofs:end-1-ofs)).^2 - tf1(:,1+ofs:end-ofs);
        title_tf = 'Time-frequency Amplitude (Taeger-Kaiser normalisation)';
    else
        power(:,2+ofs:end-1-ofs) = 0.5*abs(tf(:,2+ofs:end-1-ofs)).^2;
        title_tf = 'Time-frequency Amplitude (no normalisation)';
    end
end

function OPTIONS = select_frequency_band(fqmin, fqmax, OPTIONS)
% selection of frequency band inside the frequency bounds:
    OPTIONS.wavelet.freqs_analyzed = real(OPTIONS.wavelet.freqs_analyzed);
    Nf = length(OPTIONS.wavelet.freqs_analyzed);
    OPTIONS.wavelet.freqWindow = [1 Nf];
    a = find(OPTIONS.wavelet.freqs_analyzed<fqmin,1,'last');
    if ~isempty(a), OPTIONS.wavelet.freqWindow(1) = a; end
    a = find(OPTIONS.wavelet.freqs_analyzed>fqmax,1,'first');
    if ~isempty(a), OPTIONS.wavelet.freqWindow(2) = a; end

end


function hfig = displayTF_Plane(power,time, OPTIONS, hfig)

    if nargin  < 4
        hfig = figure('units','normalized','outerposition',[0 0 1 1]);
        ax = axes();
    else
    
        ax = gca;
    end

    set(hfig,'CurrentAxes',ax);


    % frequency cropping:
    if isfield(OPTIONS.wavelet,'freqWindow')
        power = power(OPTIONS.wavelet.freqWindow(1):OPTIONS.wavelet.freqWindow(2),:);
    else
        OPTIONS.wavelet.freqWindow = [1 size(power,1)];
    end

    % maximum power display
    disp(['Max. Value of this TF plane: ' num2str(sqrt(max(max(power))))]);
    OPTIONS.power = power;
    % bornes et marqueurs de l'axe frequence
    LFmin = log2(OPTIONS.wavelet.freqs_analyzed(OPTIONS.wavelet.freqWindow(1)));
    LFmax = log2(OPTIONS.wavelet.freqs_analyzed(OPTIONS.wavelet.freqWindow(2)));
    Ytic  = 2.^(fix(LFmin):fix(LFmax));
    temp  = num2str(Ytic','%4.3f');
    Yticl = mat2cell(temp,...
            ones(1, size(temp,1)),...
            size(temp,2));
    imagesc(time, log2(OPTIONS.wavelet.freqs_analyzed(OPTIONS.wavelet.freqWindow(1):OPTIONS.wavelet.freqWindow(2))), power);
    if isfield (OPTIONS.wavelet,'display') && isfield(OPTIONS.wavelet.display,'fontscale')
        font_scale = OPTIONS.wavelet.display.fontscale;
    else
        font_scale = 14;
    end
    
    ylabel('Frequency (Hz)')
    xlabel('Times (s)');

    set(ax,...
                'XLim',[time(1),time(end)],...
                'YGrid','on',...
                'YLim',[LFmin,LFmax],...
                'YDir','normal',...
                'YTick',log2(Ytic(:)),...
                'YTickLabel',Yticl,...
                'FontName','Times',...
                'FontAngle','Italic',...
                'FontSize',font_scale);

    if isfield(OPTIONS, 'clim') && ~isempty(OPTIONS.clim)
        clim(OPTIONS.clim);
    end

    if isfield(OPTIONS ,'colormap') && ~isempty(OPTIONS.colormap)
        colormap(OPTIONS.colormap);
    end
    grid off;
    title(ax,OPTIONS.title_tf);
end

function [spectrum_mean, spectrum_std, N] = averageBySleepStage(power,time, Events,motion)

    idx_motion      =  getIndexEvent(time,motion);
    spectrum_mean   = zeros( length(Events), size(power,1) );
    spectrum_std    = zeros( length(Events), size(power,1) );
    N               = zeros(1,length(Events));

    for iEvent = 1:length(Events)
        idx_event  =  ~idx_motion & getIndexEvent(time,Events(iEvent));
        if any(idx_event)
            spectrum_mean(iEvent,:) = mean(power(:,idx_event),2);
            spectrum_std(iEvent,:) = std(power(:,idx_event),[],2);

            N(iEvent) = sum(idx_event); 
        end
    end
end

function segments = exctractSegment(WData,time, stages, events, freqs)
% Extract continous segment of data with the same sleep stage. 
% Output :  1xN stuct segment. N being the number of segments
%      segments(i).label: stage present in the segment
%      segment(i).WDdata : tine-frequency data
%      segment(i).time   ; time starting at 0 for each segment
%      segment(i).freq   : frequency axis
%      segment(i).offset : time between the start of the recording and the
%      begining of the segment
%      segment(i).nAvg ; number of average done accros each dimension:
%      [space, time within run, time between runs  ] 
%      segment(i).duration : duration of the segment
%      segment(i).events: events occuring during the segment



% create hypnogram 
    hypnoram = zeros(1, length(time));
    for iStage = 1:length(stages)
        idx_event  =  getIndexEvent(time,stages(iStage));
        hypnoram(idx_event) = iStage;
    end
    
    % detect changes indicating begening / end of segment 
    
    diff = [ 1 , hypnoram(2:end-1) ~= hypnoram(1:end-2) , 1 ];
    changepoints = find(diff);
    
    nSegment = length(changepoints) - 1;
    segments = repmat(struct('label', '', 'WData', [], 'time' , [],'freq',[], 'offset', 0 ,'nAvg', nan(1,4), 'duration', 0 , 'events',repmat(db_template('event'),1,0) ) ,  ...
                      nSegment , 1);
    
    good_segments = true(1, nSegment );
    for iSegment = 1:nSegment
        beging = changepoints(iSegment);
        finish = changepoints(iSegment + 1) - 1;
    
        iStage = unique(hypnoram(beging:finish));


        assert(length(iStage) == 1, 'Something wrong happened, multiple sleep stage in the segment');

        if iStage == 0 % we dont have reliable info
            good_segments(iSegment) = false;
            continue;
        end

        segments(iSegment).label    = stages(iStage).label;
        segments(iSegment).WData    = WData(:, beging:finish);
        segments(iSegment).time     = time(beging:finish) - time(beging);
        segments(iSegment).freq     = freqs;
        segments(iSegment).offset   = time(beging);
        segments(iSegment).duration = time(finish) - time(beging); 
        segments(iSegment).events   = events;

        
        % we remove segment less than 2 sample
        if length(segments(iSegment).time) < 2
            good_segments(iSegment) = false;
            continue;
        end

        % fix the events : keep only the events happening in the segment
        % and adjust the time 

        good_event = true(1, length(events));

        for iEvent = 1:length(events)
            sEvent = events(iEvent);
            sEvent.times =  sEvent.times - segments(iSegment).offset;

            if isempty(sEvent.times)
                good_event(iEvent) = false;
                continue;
            end

            if size(sEvent.times, 1 )   == 1 % single events 
                inside = sEvent.times > 0 & sEvent.times < segments(iSegment).duration;
            else
                beging_inside = sEvent.times(1,:) > segments(iSegment).time(2) & sEvent.times(1,:)  < segments(iSegment).time(end-1);
                end_inside = sEvent.times(2,:) > segments(iSegment).time(2) & sEvent.times(2,:)  < segments(iSegment).time(end-1);

                inside = beging_inside | end_inside ;
                
                % we then fix the events that are at the interface
                to_fix = (beging_inside & ~end_inside);
                sEvent.times(2,to_fix) = segments(iSegment).duration;

                to_fix = (~beging_inside & end_inside);
                sEvent.times(1,to_fix) = 0;
            end

            if ~any(inside)
                good_event(iEvent) = false;
                continue;
            end

             sEvent.times  =  sEvent.times(:, inside);
             segments(iSegment).events(iEvent) = sEvent;
        end
        segments(iSegment).events = segments(iSegment).events(good_event);

    end
    
    segments = segments(good_segments);
end

function segments = epochSegment(segment, epochLength, epochOverlap)
% epochSegment: epoch the segments into sub-segment of equal lengh based on
% epochLenght. Consecutive sub-segment will overlap based on epochOverlap.
% both expressed in second. Inspired on bst_epoching
% Todo: Simplify with exctractSegment. 

    if length(segment) > 1
        segments = [];
        for kSegment = 1:length(segment)
            segments = [segments; ... 
                            epochSegment(segment(kSegment), epochLength, epochOverlap)];
        end

        return
    end
    % init variables
    fs = round(1/(segment.time(2) - segment.time(1)));
    nSamples = length(segment.time);
    % convert to number of sample
    epochLength = fs*epochLength + 1;
    epochOverlap = fs*epochOverlap;
    % Epoch shift
    nShift = (epochLength - epochOverlap);
    % Number of epochs
    nSegment = floor( (nSamples - epochOverlap) / (epochLength - epochOverlap) );
    
    if nSegment == 0
        segments = [];
        return
    end

    % Epoch Start and End indices
    epochIxs = zeros(nSegment, 2);
    epochIxs(:, 1) = ((0 : (nSegment-1))' * nShift) + 1;
    epochIxs(:, 2) = epochIxs(:, 1) + epochLength - 1;
    
    segments = repmat(segment, nSegment , 1);
    good_segments = true(1, nSegment );
    for iSegment = 1:nSegment

        beging = epochIxs(iSegment,1);
        finish = epochIxs(iSegment,2);

        segments(iSegment).WData    = segment.WData(:, beging:finish);
        segments(iSegment).time     = segment.time(beging:finish) - segment.time(beging);
        segments(iSegment).duration = round(segments(iSegment).time(end),2);
        segments(iSegment).offset   = segments(iSegment).offset + segment.time(beging);
        segments(iSegment).events   = []; % todo
    end

end


function averaged_segments = averageWithinSegment(segments)
    % Average accross time within each segment

    nSegment = length(segments);
    averaged_segments = segments;
    for iSegment = 1:nSegment
         averaged_segments(iSegment).WData = squeeze(mean(averaged_segments(iSegment).WData, 2));
         averaged_segments(iSegment).WDataStd = squeeze(std(averaged_segments(iSegment).WData,[], 2));
         averaged_segments(iSegment).nAvg(2) = length(averaged_segments(iSegment).time);
    end

    % remove time information 
    averaged_segments = rmfield(averaged_segments,"time");
    averaged_segments = rmfield(averaged_segments,"events");
    averaged_segments = rmfield(averaged_segments,"offset");

end

function resampled_segments = resampleFrequency(segments, new_f)
% Resample the frequency axis to standardized frequency axis. Mandatory to
% average between run / subjects. Only works removing the time axis for
% now (ie. after averageWithinSegment)

    nSegment = length(segments);
    resampled_segments = segments;
    for iSegment = 1:nSegment
        nTime = size(resampled_segments(iSegment).WData,2);
        resampled_segments(iSegment).WData = zeros(length(new_f), nTime);
        if isfield(resampled_segments,'WDataStd')
            resampled_segments(iSegment).WDataStd = zeros(length(new_f), nTime);
        end
        
        for iTime = 1:nTime
            resampled_segments(iSegment).WData(:,iTime) = interp1( segments(iSegment).freq, segments(iSegment).WData(:,iTime), new_f,"linear");
            if isfield(resampled_segments,'WDataStd')
                resampled_segments(iSegment).WDataStd(:,iTime) = interp1( segments(iSegment).freq, segments(iSegment).WDataStd(:,iTime), new_f,"linear");
            end
        end
        resampled_segments(iSegment).freq  = new_f;
    end

end


function averaged_segments = averageBetweenSegment(segments)

    averaged_segments = segments(1);

    wData = cat(3,segments.WData);
    averaged_segments.WData = squeeze(mean(wData,3)) ; 

    if isfield(segments(1),'WDataStd')
        disp('Replacing WDataStd by the std between segment. ')
    end
    averaged_segments.WDataStd = squeeze(std(wData,[],3));

    % Update the number of average. 

    nAvg = vertcat(segments.nAvg);
    averaged_segments.nAvg = mean(nAvg);
    averaged_segments.nAvg(3) = length(segments);

    % Update the number of average. 
    
    averaged_segments.duration = mean([segments.duration]);
end

function segments = setTimeOrigin(segments, new_zero)
    segments.time = segments.time - new_zero;
    for iEvent = 1:length(segments.events)
        segments.events(iEvent).times = segments.events(iEvent).times - new_zero;
    end
end



function f = displayPowerSpectrum(spectrum_mean,spectrum_err, labels, freqs_analyzed, OPTIONS)

    f = figure('units','normalized','outerposition',[0 0 1 1]);
    ax = axes();
    set(f,'CurrentAxes',ax);

     % bornes et marqueurs de l'axe frequence
    LFmin = log2(freqs_analyzed(1));
    LFmax = log2(freqs_analyzed(end));

    Ytic  = 2.^(fix(LFmin):fix(LFmax));
    temp  = num2str(Ytic','%4.3f');
    Yticl = mat2cell(temp,...
            ones(1, size(temp,1)),...
            size(temp,2));
    if ~isfield(OPTIONS, 'color_map')
        OPTIONS.color_map = summer(size(spectrum_mean,1));
    end

    for i = 1:size(spectrum_mean,1)
        shadedErrorBar(freqs_analyzed, spectrum_mean(i,:),spectrum_err(i,:) ,'lineProps',{'LineWidth',2,  'color',OPTIONS.color_map(i,:)}); hold on;
    end

    h = legend(ax,labels);
    title(h, 'Sleep Stage')
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax, 'Power');
    title(ax,OPTIONS.title_tf);

    set(ax,...
            'LineWidth', 2, ...
            'xscale', 'log' , ...
            'xlim', [ min(freqs_analyzed), max(freqs_analyzed)], ...
            'XTick',Ytic(:),...
            'XTickLabel',Yticl,...
            'FontName','Times',...
            'FontAngle','Italic',...
            'FontSize',  OPTIONS.wavelet.display.fontscale);

    if strcmp(OPTIONS.wavelet.display.TaegerK ,'no')
        set(ax,'yscale','log');
    end

end

