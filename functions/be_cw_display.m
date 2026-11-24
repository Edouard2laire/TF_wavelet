function hFig = be_cw_display(OPTIONS, fs, duration, fDisplay)
% be_cw_display: Interactive display of the time-frequency properties of the
% Morse wavelet filter bank built by be_cw_init.
%
%   INPUTS:
%       - OPTIONS  : options structure (same as for be_cw_init). The field
%                    OPTIONS.mandatory.DataTime is not needed, it is built
%                    here from fs and duration.
%       - fs       : sampling rate in Hz
%       - duration : (optional) file duration in seconds. If missing or empty,
%                    the user is asked for it in a dialog box.
%       - fDisplay : (optional) initial frequency in Hz at which the wavelet
%                    is drawn. Default: geometric mean of the analyzed frequencies.
%
%   OUTPUT:
%       - hFig     : handle of the figure ([] if the user cancelled)
%
% (rest of the header unchanged)

    if nargin < 3, duration = []; end
    if nargin < 4, fDisplay = []; end
    hFig = [];

    % --- Ask for the file duration if it was not provided
    if isempty(duration)
        answer = inputdlg('File duration (seconds):', 'Morse wavelet', 1, {'600'});
        if isempty(answer)
            return;                                  % cancelled
        end
        duration = str2double(answer{1});
    end
    nTime = round(duration * fs);
    if ~isscalar(duration) || ~isfinite(duration) || duration <= 0 || nTime < 8
        error('be_cw_display: invalid file duration (%s s) for a sampling rate of %g Hz.', ...
              num2str(duration), fs);
    end

    % --- Time vector consistent with the data that will be transformed
    OPTIONS.mandatory.DataTime = (0:nTime-1) / fs;

    % --- Build the filter bank (scales ordered as OPTIONS.wavelet.scale_analyzed)
    [Psi, OPTIONS] = be_cw_init(nTime, OPTIONS);
    P  = squeeze(gather(Psi));                    % nTime x (J+1)
    f  = OPTIONS.wavelet.freqs_analyzed(:);       % (J+1) x 1, in Hz, ascending
    dt = diff(OPTIONS.mandatory.DataTime([1 2]));
    N  = nTime;
    K  = numel(f);

    % --- Axes
    nHalf = floor(N/2) + 1;
    fAxis = (0:nHalf-1)' / (N*dt);                % Hz, positive frequencies
    t     = dt * (-floor(N/2) : ceil(N/2)-1)';    % s, zero lag at the center

    % --- Wavelets in the time domain (analytic, so abs() is the envelope)
    W   = fftshift(ifft(P, [], 1), 1);            % nTime x (J+1)
    env = abs(W);
    amp = abs(P(1:nHalf, :));

    % --- Resolutions (FWHM of the amplitude)
    FWHM_f = zeros(K,1);
    FWHM_t = zeros(K,1);
    for k = 1:K
        FWHM_f(k) = compute_fwhm(fAxis, amp(:,k));
        FWHM_t(k) = compute_fwhm(t,     env(:,k));
    end

    % --- Initial wavelet
    if isempty(fDisplay)
        fDisplay = sqrt(min(f) * max(f));
    end
    [~, iDisp] = min(abs(f - fDisplay));

    % --- Figure
    hFig = figure(...
        'MenuBar',     'none', ...
        'Toolbar',     'none', ...
        'NumberTitle', 'off', ...
        'Name',        'Morse wavelet', ...
        'Pointer',     'arrow');
    hAxes1 = axes();
    hAxes2 = axes();
    hAxes3 = axes();
    fontSize = bst_get('FigFont');

    % Spectral resolution + marker
    plot(hAxes1, f, FWHM_f);
    hold(hAxes1, 'on');
    hMark1 = plot(hAxes1, f(iDisp), FWHM_f(iDisp), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    title(hAxes1, 'Spectral resolution', 'Interpreter', 'none', 'FontSize', fontSize);
    xlabel(hAxes1, 'Frequency (Hz)', 'FontSize', fontSize);
    ylabel(hAxes1, 'FWHM (Hz)', 'FontSize', fontSize);
    set(hAxes1, 'Position', [.08, .64, .34, .28]);

    % Temporal resolution (log frequency axis) + marker
    semilogx(hAxes2, f, FWHM_t);
    hold(hAxes2, 'on');
    hMark2 = plot(hAxes2, f(iDisp), FWHM_t(iDisp), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 7);
    title(hAxes2, 'Temporal resolution', 'Interpreter', 'none', 'FontSize', fontSize);
    xlabel(hAxes2, 'Frequency (Hz)', 'FontSize', fontSize);
    ylabel(hAxes2, 'FWHM (sec)', 'FontSize', fontSize);
    set(hAxes2, 'Position', [.08, .28, .34, .28], 'XLim', [min(f), max(f)]);

    % Wavelet (graphics objects are created once, then updated by the slider)
    hReal = plot(hAxes3, t, real(W(:,iDisp)), 'linewidth', 2);
    hold(hAxes3, 'on');
    hImag = plot(hAxes3, t, imag(W(:,iDisp)), 'r', 'linewidth', 2);
    hLineLo = line([0 0], [-1 1], 'LineWidth', 1, 'Color', [.7 .7 .7], 'Parent', hAxes3);
    hLineHi = line([0 0], [-1 1], 'LineWidth', 1, 'Color', [.7 .7 .7], 'Parent', hAxes3);
    hText = text(0, 1, '99% energy', ...
        'Color',               [.5 .5 .5], ...
        'FontSize',            fontSize, ...
        'FontUnits',           'points', ...
        'VerticalAlignment',   'top', ...
        'HorizontalAlignment', 'right', ...
        'Parent',              hAxes3);
    hTitle = title(hAxes3, '', 'Interpreter', 'none');
    xlabel(hAxes3, 'Time (s)', 'FontSize', fontSize);
    ylabel(hAxes3, 'Amplitude', 'FontSize', fontSize);
    set(hAxes3, 'Position', [.50, .30, .44, .62]);

    % Fonts and grids
    set([hAxes1, hAxes2, hAxes3], ...
        'FontSize', fontSize, ...
        'XGrid', 'on', 'YGrid', 'on', 'XMinorGrid', 'on', 'YMinorGrid', 'on', ...
        'XMinorTick', 'on', 'YMinorTick', 'on');

    % --- Slider (one step per analyzed frequency, i.e. log-spaced in frequency)
    hSlider = uicontrol(hFig, 'Style', 'slider', 'Units', 'normalized', ...
        'Position', [.50, .20, .30, .035], ...
        'Min', 1, 'Max', K, 'Value', iDisp, ...
        'SliderStep', [1/(K-1), min(1, 5/(K-1))]);
    hLabel = uicontrol(hFig, 'Style', 'text', 'Units', 'normalized', ...
        'Position', [.81, .19, .14, .05], 'FontSize', fontSize, ...
        'HorizontalAlignment', 'left');

    % Description
    strDesc = ['The complex Morse wavelet is analytic: its spectrum is zero at negative ' ...
               'frequencies (blue: real part, red: imaginary part). Its shape is not Gaussian, ' ...
               'so resolutions are measured numerically as the FWHM (full width half maximum) ' ...
               'of the amplitude spectrum and of the temporal envelope, at each analyzed frequency. ' ...
               'Use the slider to move the red markers and display the wavelet at another frequency.'];
    annotation(hFig, 'textbox', [.05, .02, .90, .10], 'String', strDesc, 'EdgeColor', 'none');

    % Continuous update while dragging (fallback to the callback on older releases)
    addlistener(hSlider, 'ContinuousValueChange', @(~,~) updateDisplay());
    set(hSlider, 'Callback', @(~,~) updateDisplay());

    updateDisplay();

    % ====================================================================
    function updateDisplay()
        k = min(max(round(get(hSlider, 'Value')), 1), K);
        w = W(:, k);

        % Markers on the resolution plots
        set(hMark1, 'XData', f(k), 'YData', FWHM_f(k));
        set(hMark2, 'XData', f(k), 'YData', FWHM_t(k));

        % Interval containing the central 99% of the wavelet energy
        cumE = cumsum(abs(w).^2) / sum(abs(w).^2);
        tLo  = t(find(cumE >= 0.005, 1, 'first'));
        tHi  = t(find(cumE >= 0.995, 1, 'first'));
        tLim = 1.5 * max(abs([tLo, tHi]));
        if isempty(tLim) || tLim == 0
            tLim = max(abs(t));
        end
        yLim = [-1, 1] * 1.05 * max(abs(w));

        % Wavelet
        set(hReal, 'YData', real(w));
        set(hImag, 'YData', imag(w));
        set(hAxes3, 'XLim', [-tLim, tLim], 'YLim', yLim);
        set(hLineLo, 'XData', tLo*[1 1], 'YData', yLim);
        set(hLineHi, 'XData', tHi*[1 1], 'YData', yLim);
        set(hText, 'Position', [tHi - 0.02*tLim, yLim(2), 0]);
        set(hTitle, 'String', sprintf('Complex Morse wavelet (n=%g, m=%g) at %.3g Hz', ...
            OPTIONS.wavelet.vanish_moments, OPTIONS.wavelet.order, f(k)));

        % Label
        set(hLabel, 'String', sprintf('%.3g Hz (period %.3g s)', f(k), 1/f(k)));
    end
end

% ========================================================================
function w = compute_fwhm(x, y)
% Full width at half maximum of a single-peaked curve y(x), using linear
% interpolation at the two half-maximum crossings. NaN if a crossing is not
% reached inside the axis.
    [ymax, ip] = max(y);
    h  = ymax / 2;
    iL = find(y(1:ip) < h, 1, 'last');
    iR = ip - 1 + find(y(ip:end) < h, 1, 'first');
    if isempty(iL) || isempty(iR) || iL >= ip || iR <= ip
        w = NaN;
        return;
    end
    xL = x(iL)   + (h - y(iL))   / (y(iL+1) - y(iL))   * (x(iL+1) - x(iL));
    xR = x(iR-1) + (h - y(iR-1)) / (y(iR)   - y(iR-1)) * (x(iR)   - x(iR-1));
    w  = xR - xL;
end