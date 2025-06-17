function [power, title_tf] = be_apply_measure(tf, OPTIONS)

    ofs     = ceil(0.02*size(tf,3));

    power   = zeros(size(tf));
    power(:,:,2+ofs:end-1-ofs) = 0.5*abs(tf(:,:,2+ofs:end-1-ofs)).^2;

    title_tf = 'Time-frequency Amplitude';

    if isfield(OPTIONS.wavelet.display,'TaegerK') && strcmp(OPTIONS.wavelet.display.TaegerK,'yes') 

        for t = (2+ofs):( size(tf,3)-1-ofs)
            tf1 = 0.25*tf(:,:,t+1).*conj(tf(:,:,t-1))+0.25*conj(tf(:,:,t+1)).*tf(:,:,t-1);       
            power(:,:,t) = power(:,:,t)  - tf1;
        end

        title_tf = [title_tf, '(Taeger-Kaiser normalisation)'];
    end

end