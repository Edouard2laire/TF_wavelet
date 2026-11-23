function [power, title_tf] = be_apply_measure(tf, OPTIONS)

    ofs     = ceil(0.02*size(tf,3));

    power   = zeros(size(tf));
    power(:,:,2+ofs:end-1-ofs) = 0.5*abs(tf(:,:,2+ofs:end-1-ofs)).^2;

    title_tf = 'Time-frequency Amplitude';

    if isfield(OPTIONS.wavelet.display,'TaegerK') && strcmp(OPTIONS.wavelet.display.TaegerK,'yes') 

        % Teager-Kaiser energy operator: 0.5*( |x(t)|^2 - Re[ x(t+1) * conj(x(t-1)) ] )
        iTime      = (2+ofs):(size(tf,3)-1-ofs);
        teagerTerm = 0.5 * real( tf(:,:,iTime+1) .* conj(tf(:,:,iTime-1)) );
        power(:,:,iTime) = power(:,:,iTime) - teagerTerm;

        title_tf = [title_tf, '(Taeger-Kaiser normalisation)'];
    end

end