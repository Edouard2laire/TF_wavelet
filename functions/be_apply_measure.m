function [power, title_tf] = be_apply_measure(tf, OPTIONS)

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