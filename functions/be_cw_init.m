function [Psi, OPTIONS] = be_cw_init(nTime, OPTIONS)
% be_cw_init: Precompute the Morse wavelet filter bank in the Fourier domain
% and the wavelet analysis options. Call it once, then reuse Psi for every
% call to be_cw_transform.
%
%   INPUTS:
%       - nTime   : number of time samples of the data that will be transformed
%       - OPTIONS : options structure
%
%   OUTPUTS:
%       - Psi     : real array, size [1 x nTime x (J+1)], scales ordered as
%                   OPTIONS.wavelet.scale_analyzed (from largest to smallest scale)
%       - OPTIONS : options structure, with the fields scale_analyzed,
%                   freqs_analyzed and Scone filled in
%
% (Copyright / license header: same as be_cw_transform)

    % Set wavelet parameters:
    n  = OPTIONS.wavelet.vanish_moments;
    m  = OPTIONS.wavelet.order;
    J  = OPTIONS.wavelet.nb_levels;
    dt = diff( OPTIONS.mandatory.DataTime([1 2]) );
    N  = nTime;

    % --- definition de l'ondelette
    [wmin, wmax] = Morse_support_spectral(n,m);
    [tmin, tmax] = Morse_support_temporel(n,m);

    % frequence centrale de l'ondelette de Cauchy
    f0 = (((n+0.5)/m).^(1/m))/2/pi;

    % --- echelles min et max
    a_max = (N-3)*dt/(tmax-tmin);
    a_min = wmax/pi*dt;

    % --- resolution en octaves et en voies
    dj = log2(a_max/a_min)/J;
    aj = a_min*2.^((0:J)*dj);

    % --- options
    OPTIONS.wavelet.scale_analyzed = fliplr(aj);
    OPTIONS.wavelet.freqs_analyzed = f0./fliplr(aj);
    OPTIONS.wavelet.Scone          = (f0/dt)./[1E-9,1:((N+1)/2-1),fliplr((1:(N/2-1))),1E-9];

    % --- espace de Fourier
    wk = 2.*pi/N/dt*(1:fix(N/2));
    wk = [0., wk, -wk(fix((N-1)/2):-1:1)];

    % --- filter bank: (J+1) x N, then 1 x N x (J+1)
    Psi = Morse_w(aj' * wk, n, m) .* sqrt(aj');
    Psi = permute(Psi, [3 2 1]);

    % Flip the scale dimension once here, so be_cw_transform does not have to
    Psi = flip(Psi, 3);
end

% ========================================================================
function y = Morse_w(w,n,m)
    y  = zeros(size(w));
    wp = w(w>0);
    y(w>0) = wp.^n.*exp(-wp.^m);
end

function [tmin,tmax] = Morse_support_temporel(n,m)
    s = 0.01;
    tmax = sqrt(s^(-1/(n+1))*(gamma(n+1)/2/pi)^(2/(n+1))-1);
    tmin = -tmax;
end

function [wmin,wmax] = Morse_support_spectral(n,m)
    s  = 0.01;
    s0 = s^(1/n);
    s1 = -log(s)/n;
    w0 = 0;
    for i=1:50
        w0 = w0^m;
        w0 = s0*exp(w0/n);
    end
    wmin = w0;
    w0 = n;
    for i = 1:50
        w0 = s1+n/m*log(w0);
    end
    wmax = w0^(1/m);
end