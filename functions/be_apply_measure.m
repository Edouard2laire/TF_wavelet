function [power, title_tf] = be_apply_measure(tf, OPTIONS)
% be_apply_measure: Compute the time-frequency power map from wavelet coefficients.
%
% Inputs:
%   tf      - Complex time-frequency (wavelet) coefficients, size [nChan x nFreq x nTime]
%   OPTIONS - Options structure; OPTIONS.wavelet.measure selects the
%             energy computation ['none', 'power', 'TeagerK'] (see below)
%
% Output:
%   power   - Real array, same size as tf. The first and last (1+ofs) time
%             samples are left at zero.
%
% Three modes, selected by OPTIONS.wavelet.display.TaegerK:
%   0) 'none', do nothing. Return the complex values
%
%   1) 'power': squared magnitude of the coefficients
%        power(t) = 0.5 * |x(t)|^2
%
%   2) 'TeagerK' (default): Teager-Kaiser energy operator (TKEO) along the time axis,
%        power(t) = 0.5 * |x(t)|^2 - 0.5 * Re[ x(t+1) * conj(x(t-1)) ]
%      This tracks the instantaneous energy of the signal (amplitude and frequency) 
%      and is more sensitive to short transients than |x(t)|^2; notably
%      removing the effect of the 1/f in the time-frequency plane
%
% where x(t) = tf(:,:,t) and t spans (2+ofs):(end-1-ofs).
    
    % Initialize output
    measure   = OPTIONS.wavelet.measure;
    ofs       = ceil(0.02*size(tf,3));
    iTime     = (2+ofs):(size(tf,3)-1-ofs);
    power     = zeros(size(tf));

    % If measure if none, do nothing
    if strcmpi(measure, 'none')
        power(:, :, iTime) = tf(:, :, iTime);
        title_tf = '';
        return
    end 

    useTeager = strcmpi(measure, 'TeagerK');
    if useTeager
        % power(t) = 0.5 * |x(t)|^2 - 0.5 * Re[ x(t+1) * conj(x(t-1)) ]
        power(:,:,iTime) = 0.5*abs(tf(:,:,iTime)).^2 - 0.5 * real(tf(:,:,iTime+1) .* conj(tf(:,:,iTime-1)));
        title_tf = 'Time-frequency Power (Teager-Kaiser energy operator)';
    else
        % power(t) = 0.5 * |x(t)|^2 
        power(:,:,iTime) = 0.5*abs(tf(:,:,iTime)).^2;
        title_tf = 'Time-frequency Power';
    end
end