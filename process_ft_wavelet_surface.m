
function varargout = process_ft_wavelet_surface( varargin )
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
% Authors: Jean-Marc Lina, Edouard Delaire (2023)

eval(macro_method);
end


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription()
    % Description the process
    sProcess.Comment     = 'Time-frequency: Morse Wavelet(source)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = {'NIRS','Frequency'};
    sProcess.Index       = 1705;
    sProcess.Description = 'http://www.fieldtriptoolbox.org/tutorial/timefrequencyanalysis';
    % Definition of the input accepted by this process
    sProcess.InputTypes  = {'import', 'data'};
    sProcess.OutputTypes = {'import', 'data'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Definition of the input accepted by this process
    sProcess.InputTypes  = { 'results'};
    sProcess.OutputTypes = { 'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % Options: Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    sProcess.options.timewindow.Group   = 'input';
    

    % === CLUSTERS
    sProcess.options.scouts.Comment = '';
    sProcess.options.scouts.Type    = 'scout';
    sProcess.options.scouts.Value   = {};


    sProcess.options.label1.Comment = '<b>Wavelet parameters:</b>';
    sProcess.options.label1.Type    = 'label';


    sProcess.options.vanish_moments.Comment = 'Vanishing Moment: ';
    sProcess.options.vanish_moments.Type    = 'value';
    sProcess.options.vanish_moments.Value   = {1.5, '',1};

    sProcess.options.order.Comment = 'Order: ';
    sProcess.options.order.Type    = 'value';
    sProcess.options.order.Value   = {20, '',0};

    sProcess.options.nb_levels.Comment = 'Nb levels: ';
    sProcess.options.nb_levels.Type    = 'value';
    sProcess.options.nb_levels.Value   = {1024, '',0};

    sProcess.options.label2.Comment = '<b>Display options:</b>';
    sProcess.options.label2.Type    = 'label';


    sProcess.options.normalization.Comment = 'Apply normalization: ';
    sProcess.options.normalization.Type    = 'checkbox';
    sProcess.options.normalization.Value   = 1;


end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess)
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs)
    OutputFiles = {};
    
    bst_plugin('Load', 'tfnirs');

    % Load ROIs
    sSubject    = bst_get('Subject',sInputs(1).SubjectName );
    sCortex     = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);

    % Load Rois
    [tempalte_OPTIONS, selected_scouts] = process_ft_wavelet('getOptions', sProcess);
    selected_scouts = getScouts(sCortex, selected_scouts{1}, selected_scouts{2});
    
    % Get Memory information
    ProcessOptions = bst_get('ProcessOptions');
    MaxSizeDouble = ProcessOptions.MaxBlockSize;

    for iCluster = 1:length(selected_scouts) 

        sROI    = selected_scouts(iCluster);
        vertex  = sROI.Vertices;
        
        ref = [];

        for iFile = 1:length(sInputs)

            OPTIONS = tempalte_OPTIONS;
            sData   = in_bst_results(sInputs(iFile).FileName, 0);
 
            sizePerSample    = prod([length(sData.Time), 1+OPTIONS.wavelet.nb_levels]);
            nSamples = length(vertex);

            blockSize   = min(max(floor(MaxSizeDouble / sizePerSample), 1), nSamples);
            nBlocks     = ceil(nSamples / blockSize);

            bst_progress('start', 'Running Time-Frequency Analysis', 'Running Time-Frequency Analysis', 0, nBlocks);

            [~, iOrigTime, ~]  = getData(sData , OPTIONS, vertex(1));
            meanvar_estimator  = be_WelfordVariance(nSamples,  [1+OPTIONS.wavelet.nb_levels, sum(iOrigTime)]);
    
            for iBlock = 1 : nBlocks

                idx = [((iBlock-1) * blockSize + 1) : min(iBlock * blockSize, nSamples)];
                [data_cortex, iOrigTime, OPTIONS]  = getData(sData , OPTIONS, vertex(idx));
    
                % Step 1 - compute time-frequency representation
                fprintf('Tf-nirs [%d/%d] > Computing the CW decomposition ... ', iBlock, nBlocks)
                [power, OPTIONS] = be_CWavelet(data_cortex, OPTIONS);
                power = permute(power(:, iOrigTime, :), [1, 3, 2]); % nVertex x Nfreq x Ntime
                fprintf(' Done. \n')
    
                 % Step 2- Normalize the TF maps
                fprintf('Tf-nirs [%d/%d]> Computing the Teager Keaser nornmalization ... ', iBlock, nBlocks)
                [power, title_tf] = be_apply_measure(power, OPTIONS);
                fprintf(' Done. \n')
        
                % Step 3 - Average accross vertex
                meanvar_estimator.update(power);
                bst_progress('inc', 1);
            end
            wData = reshape(meanvar_estimator.Mean, [1, 1+OPTIONS.wavelet.nb_levels, sum(iOrigTime)]);
            OPTIONS.title_tf =  title_tf ;

            % Step 4- Normalize the TF maps ( standardize power)
            fprintf('Tf-nirs > Normalizing the power of each vertex... ')
            [wData, ref] = be_scale_TF(wData, ref);
            wData = squeeze(wData);
            fprintf(' Done. \n')
    
            OPTIONS.title_tf  = sprintf('%s - %s', strrep(sROI.Label,'_',' '), OPTIONS.title_tf);
    
            nfreq = length(OPTIONS.wavelet.freqs_analyzed);
            nTime = length(sData.Time);
            
            TFmask  = ones(nfreq, nTime);
            is_zero = all(squeeze(power(1,:, :)) == 0);
            TFmask(:, is_zero) = 0;
    
            sDataOut            = db_template('timefreqmat');
            sDataOut.TF(1, :, :) = wData';
            sDataOut.RowNames   =  {sROI.Label};
            sDataOut.TFmask     = TFmask;
            sDataOut.DataType   = 'scout';
            sDataOut.DataFile   = sInputs(iFile).FileName;
            sDataOut.Time       = sData.Time;
            sDataOut.Options    = OPTIONS;
            sDataOut.Freqs      = OPTIONS.wavelet.freqs_analyzed;
            sDataOut.Comment    = sprintf('CW [%s]', sROI.Label);
            sDataOut.Method     = 'morlet';
            sDataOut.Measure    = 'power';
            sDataOut.Options.PowerUnits = 'physical';
            sDataOut.DisplayUnits       =' \mumol.l-1'; 
    
            % Generate a new file name in the same folder
            sStudy = bst_get('Study', sInputs(iFile).iStudy);
            OutputFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'timefreq_wavelet');
            sDataOut.FileName = file_short(OutputFile);
            bst_save(OutputFile, sDataOut, 'v7');
            % Register in database
            db_add_data(sInputs(iFile).iStudy, OutputFile, sDataOut);
            OutputFiles{end+1} = OutputFile;
        end
    end

    bst_progress('stop');
end 

function selected_scouts = getScouts(sCortex, selected_atlas, selected_scouts)

    iAtlas      = find(strcmp( {sCortex.Atlas.Name}, selected_atlas));
    iRois       = cellfun(@(x)find(strcmp( {sCortex.Atlas(iAtlas).Scouts.Label},x)), selected_scouts);
    
    selected_scouts = sCortex.Atlas(iAtlas).Scouts(iRois);
end

function [data, iOrigTime, OPTIONS]  = getData(sData , OPTIONS, vertex)

    F = sData.ImageGridAmp;

    if isempty(OPTIONS.TimeWindow)
        time = sData.Time;
        iTime = 1:length(time);
    else
        iTime = panel_time('GetTimeIndices', sData.Time, TimeWindow');
        time  = sData.Time(iTime);
    end
    OPTIONS.mandatory.DataTime     = time; 

    % Find the appropriate padding for this dataset
    p           = nextpow2(length(time)); 
    padding     = [0 round((2^p - size(time,2))/2)];
    ext_time    = padarray(time, padding, NaN,'both');
    iOrigTime   = ~isnan(ext_time);

    if ~iscell(F)
        F = {F};
    end

    % create spatial mapping and reduce the spatial dimension 
    mapping = sparse(1:length(vertex), vertex, 1, length(vertex), size(F{1},1));
    F{1}  = mapping * F{1};

    % create the temporal mapping and reduce the temporal dimension 
    mapping = sparse(iTime,  1:length(iTime), 1, size(F{end},2), length(iTime));
    F{end} = F{end} * mapping;

    % Apply the mapping to the data
    tmp = F{1};
    for iDecomposition = 2 : length(F)
     tmp = tmp * F{iDecomposition};
    end

    data = full(tmp);
    data = padarray(data, padding, 0, 'both');

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
    Event.times(1,:) = Event.times(1,:) - before;
    Event.times(2,:) = Event.times(2,:) + after;

    Event.times = Event.times(:, Event.times(2,:) > Event.times(1,:));

end

function f = displayTF_Plane(power,time, OPTIONS)

    f = figure('units','normalized','outerposition',[0 0 1 1]);
    ax = axes();
    set(f,'CurrentAxes',ax);
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
    title(ax,OPTIONS.title_tf);


end


function f = displayPowerSpectrum(power,time, Events,motion, OPTIONS)


    if isfield(OPTIONS.wavelet,'freqWindow')
        power = power(OPTIONS.wavelet.freqWindow(1):OPTIONS.wavelet.freqWindow(2),:);
        freqs_analyzed = OPTIONS.wavelet.freqs_analyzed(OPTIONS.wavelet.freqWindow(1):OPTIONS.wavelet.freqWindow(2));
    else
        OPTIONS.wavelet.freqWindow = [1 size(power,1)];
        freqs_analyzed =  OPTIONS.wavelet.freqs_analyzed;
    end
    
    idx_motion      =  getIndexEvent(time,motion);
    spectrum_mean   = zeros( length(Events), size(power,1) );
    spectrum_std    = zeros( length(Events), size(power,1) );
    N               = zeros(1,length(Events));
    for iEvent = 1:length(Events)
        idx_event  =  ~idx_motion & getIndexEvent(time,Events(iEvent));
        if any(idx_event)
            spectrum_mean(iEvent,:) = mean(power(:,idx_event),2);
            spectrum_std(iEvent,:) = std(power(:,idx_event),[],2);

            % Try to estimate the number of independant sample
            N(iEvent) = ess(power(:,idx_event), true, 'bulk'); %sum(idx_event); 
        end
    end

    f = figure('units','normalized','outerposition',[0 0 1 1]);
    ax = axes();
    set(f,'CurrentAxes',ax);

    for i = 1:size(spectrum_mean,1)
        shadedErrorBar(freqs_analyzed, spectrum_mean(i,:),spectrum_std(i,:)./sqrt(N(i)) ,'lineProps',{'color',OPTIONS.color_map(i,:)}); hold on;
    end
    xlim( [ min(freqs_analyzed), max(freqs_analyzed)])
    set(ax,'xscale','log');
    
    if strcmp(OPTIONS.wavelet.display.TaegerK ,'no')
        set(ax,'yscale','log');
    end

    legend({Events.label})
    xlabel('Frequency (Hz)')
    ylabel('Power');
    title(OPTIONS.title_tf)
end

