function [power, title_tf] = be_apply_measure(tf, OPTIONS)
% be_apply_measure: Compute the time-frequency power map from wavelet coefficients.
%
% Inputs:
%   tf      - Complex time-frequency (wavelet) coefficients, size [nChan x nFreq x nTime]
%   OPTIONS - Options structure; OPTIONS.wavelet.display.TaegerK selects the
%             energy computation (see below)
%
% Output:
%   power   - Real array, same size as tf. The first and last (1+ofs) time
%             samples are left at zero.
%
% Two modes, selected by OPTIONS.wavelet.display.TaegerK:
%
%   1) Default (TaegerK absent or 'no'): squared magnitude of the coefficients
%        power(t) = 0.5 * |x(t)|^2
%
%   2) TaegerK = 'yes': Teager-Kaiser energy operator (TKEO) along the time axis,
%        power(t) = 0.5 * |x(t)|^2 - 0.5 * Re[ x(t+1) * conj(x(t-1)) ]
%      This tracks the instantaneous energy of the signal 
%      (amplitude andcfrequency) and is more sensitive to short transients than |x(t)|^2.
%
% where x(t) = tf(:,:,t) and t spans (2+ofs):(end-1-ofs).

    ofs     = ceil(0.02*size(tf,3));

    power   = zeros(size(tf));
    iTime      = (2+ofs):(size(tf,3)-1-ofs);
    power(:,:, iTime) = 0.5*abs(tf(:,:,iTime)).^2;
    title_tf = 'Time-frequency Amplitude';

    if isfield(OPTIONS.wavelet.display,'TaegerK') && strcmp(OPTIONS.wavelet.display.TaegerK,'yes') 

        % Teager-Kaiser energy operator: 0.5*( |x(t)|^2 - Re[ x(t+1) * conj(x(t-1)) ] )
        teagerTerm = 0.5 * real( tf(:,:,iTime+1) .* conj(tf(:,:,iTime-1)) );
        power(:,:,iTime) = power(:,:,iTime) - teagerTerm;

        title_tf = [title_tf, '(Taeger-Kaiser normalisation)'];
    end

end