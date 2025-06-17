
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
function sProcess = GetDescription() %#ok<DEFNU>
    % Description the process
    sProcess.Comment     = 'Time-frequency: Morse Wavelet(source)';
    sProcess.Category    = 'File';
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
    % Options: Sensor types

    sProcess.options.sensortypes.Comment = 'Sensor types or names (empty=all): ';
    sProcess.options.sensortypes.Type    = 'text';
    sProcess.options.sensortypes.Value   = 'NIRS';

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
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

    sData            = load(file_fullpath(sInputs.FileName));

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

    % Find the appropriate padding for this dataset
    p           = nextpow2(length(time)); 
    padding     = [0 round((2^p - size(time,2))/2)];
    ext_time    = padarray(time, padding, NaN,'both');
    iOrigTime   = ~isnan(ext_time);

    % Load data:
    F           = sData.ImageGridAmp;

    % Load ROIs
    sSubject    = bst_get('Subject',sInputs.SubjectName );
    sCortex     = in_tess_bst(sSubject.Surface(sSubject.iCortex).FileName);
    ROI         = sProcess.options.scouts.Value;
    iAtlas      = find(strcmp( {sCortex.Atlas.Name},ROI{1}));
    iRois       = cellfun(@(x)find(strcmp( {sCortex.Atlas(iAtlas).Scouts.Label},x)),   ROI{2});

    ROI_select = {sCortex.Atlas(iAtlas).Scouts(iRois).Label};
    if isempty(iRois)
        disp('ROI doesnt exist.');
        return;
    end

    OPTIONS.wavelet.vanish_moments =  sProcess.options.vanish_moments.Value{1} ; % vanish momemts
    OPTIONS.wavelet.order          =  sProcess.options.order.Value{1} ;  % spectral decay
    OPTIONS.wavelet.nb_levels      =  sProcess.options.nb_levels.Value{1};% number of voices
    OPTIONS.wavelet.verbose        = 0;   % verbose or not
    OPTIONS.mandatory.DataTime     = time; 

    if sProcess.options.normalization.Value
        OPTIONS.wavelet.display.TaegerK = 'yes';
    else
        OPTIONS.wavelet.display.TaegerK = 'no';
    end

    

    OPTIONS     = repmat(OPTIONS, 1,length(ROI_select) );
    nSensors    = sum(cellfun(@(x)length(x),{sCortex.Atlas(iAtlas).Scouts(iRois).Vertices}));

    bst_progress('start', 'Running Time-Frequency Analysis', 'Running Time-Frequency Analysis', 0, nSensors);
    for iCluster = 1:length(ROI_select) 
        sROI = sCortex.Atlas(iAtlas).Scouts(iRois(iCluster));
        vertex = sROI.Vertices;
        

        data_cortex = getData(F , vertex, iTime, padding);

        % Step 1 - compute time-frequency representation
        fprintf('Tf-nirs > Computing the CW decomposition ... ')
        [power, OPTIONS(iCluster)] = be_CWavelet(squeeze(data_cortex(:, 1:2^p)), OPTIONS(iCluster));
        power = permute(power, [1, 3, 2]); % N_channel x Nfreq x Ntime
        fprintf(' Done. \n')


        % Step 2- Normalize the TF maps
        fprintf('Tf-nirs > Computing the Teager Keaser nornmalization ... ')
        [power, title_tf] = apply_measure(power(:, :, iOrigTime), OPTIONS);
        OPTIONS(iCluster).title_tf =  title_tf ;
        fprintf(' Done. \n')
    

        % Step 3- Normalize the TF maps ( standardize power)
        fprintf('Tf-nirs > Normalizing the power of each vertex... ')
        for iSensor = 1:length(vertex)

            power_time         = sqrt(sum(power(iSensor,:,:).^2));
            power(iSensor,:,:) = power(iSensor,:,:) ./ median(power_time) ;

            bst_progress('inc', 1); 
        end
        fprintf(' Done. \n')

        % Step 3 - Average accross 
        if length(vertex) > 1
            wData =  mean(power);
        else
            wData =  power;
        end

        OPTIONS(iCluster).title_tf  = sprintf('%s - %s', strrep(sROI.Label,'_',' '), OPTIONS(iCluster).title_tf);

        nfreq = length(OPTIONS(iCluster).wavelet.freqs_analyzed);
        nTime = length(time);
        
        TFmask  = ones(nfreq,nTime);
        is_zero = all(squeeze(wData(1,:,:)) == 0);
        TFmask(:, is_zero) = 0;

        sDataOut            = db_template('timefreqmat');
        sDataOut.TF         = permute(wData, [1,3,2]); %nSensor x nTime x nFreq
        sDataOut.RowNames   =  {sROI.Label};
        sDataOut.TFmask     = TFmask;
        sDataOut.DataType   = 'scout';
        sDataOut.DataFile   = sInputs.FileName;
        sDataOut.Time       = time;
        sDataOut.Options    = OPTIONS(iCluster);
        sDataOut.Freqs      = OPTIONS(iCluster).wavelet.freqs_analyzed;
        sDataOut.Comment    = sprintf('CW [%s]', sROI.Label);
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


function data = getData(F , vertex, iTime, padding)

    if ~iscell(F)
        F = {F};
    end

    % create spatial mapping and reduce the spatial dimension 
    mapping = zeros(length(vertex), size(F{1},1)); 
    for iNode = 1:length(vertex)
        mapping(iNode,vertex(iNode)) = 1;
    end
    mapping  = sparse(mapping);
    F{1}  = mapping * F{1};

    % create the temporal mapping and reduce the temporal dimension 
    mapping = zeros(size(F{end},2), length(iTime)); 
    for iNode = 1:length(iTime)
        mapping(iTime(iNode),iNode) = 1;
    end
    mapping  = sparse(mapping);
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

function [power, title_tf] = normalize(tf, OPTIONS)
    [nCh, nT] = size(tf);
    ofs = ceil(0.02 * nT);
    idx_start = 2 + ofs;
    idx_end   = nT - 1 - ofs;
    
    % Initialize power matrix (can also be omitted if only storing the output slices)
    power = zeros(nCh, nT);
    
    % Define window parameters
    win_size = 200;  % adjust based on available memory
    step_size = win_size;  % can use overlap if needed
    
    % Loop over time in windows
    for t_start = idx_start:step_size:idx_end
        t_end = min(t_start + win_size - 1, idx_end);
        
        % Indices for tf1 = tf(:,3:end) * conj(tf(:,1:end-2));
        t1 = t_start - 1;
        t2 = t_start + 1;

        if t1 < 1 || t2 > nT
            continue;  % skip invalid bounds
        end
    
        % Extract only the needed slices
        tf_m1 = tf(:, t1:t_end-1);     % tf(:,1:end-2)
        tf_p1 = tf(:, t2:t_end+1);     % tf(:,3:end)
        tf_0  = tf(:, t_start:t_end);  % tf(:,2+ofs:end-1-ofs)
    
        % Compute tf1 and power for this chunk
        tf1_chunk = 0.25 * tf_p1 .* conj(tf_m1) + 0.25 * conj(tf_p1) .* tf_m1;
        power_chunk = 0.5 * abs(tf_0).^2 - tf1_chunk;
    
        % Store result
        power(:, t_start:t_end) = power_chunk;
        title_tf = 'Time-frequency Amplitude (Taeger-Kaiser normalisation)';
    end
end

function [power, title_tf] = apply_measure(tf, OPTIONS)

    ofs     = ceil(0.02*size(tf,3));
    power   = zeros(size(tf));

    if isfield(OPTIONS.wavelet.display,'TaegerK') && strcmp(OPTIONS.wavelet.display.TaegerK,'yes') 
        tf1 = 0.25*tf(:,:,3:end).*conj(tf(:,:,1:end-2))+0.25*conj(tf(:,:,3:end)).*tf(:,:,1:end-2);       
        power(:,:,2+ofs:end-1-ofs) = 0.5*abs(tf(:,:,2+ofs:end-1-ofs)).^2 - tf1(:,:,1+ofs:end-ofs);
        title_tf = 'Time-frequency Amplitude (Taeger-Kaiser normalisation)';
    else
        power(:,2+ofs:end-1-ofs) = 0.5*abs(tf(:,2+ofs:end-1-ofs)).^2;
        title_tf = 'Time-frequency Amplitude (no normalisation)';
    end


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

