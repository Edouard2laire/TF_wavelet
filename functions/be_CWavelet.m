function [WData, OPTIONS] = be_CWavelet(data, OPTIONS, varargin)
% cette fonction calcule la transformee en ondelette complexe du signal
% (data.Signal) echantillonne (1/data.Sampling). 
% La famille d'ondelette est celle de Morse (dont Lusin et Cauchy sont des
% cas particuliers), d'ordre n (n=1,2,...). 
% J est le nombre de voix entre l'echelle min 
% (qui doit etre superieure a l'echelle minimale ateignable, cf. a_min)
% et l'echelle max (inferieure a l'echelle maximale ateignable, cf. a_max)
%
%   INPUTS:
%       -   data    :   data matrix
%       -   OPTIONS :   options structure
%
%   OUTPUTS:
%       -   WData   :   wavelet transform
%       -   OPTIONS :   options structure
%
%% ==============================================   
% Copyright (C) 2012 - LATIS Team
%
%  Authors: JM Lina, 2012, jan 1st
%
%% ==============================================
% License 
%
% BEst is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    BEst is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with BEst. If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------   

    options_only = 0;
    if nargin==1
        options_only = 1;
    end
    
    % Set wavelet parameters:
    n       = OPTIONS.wavelet.vanish_moments;
    m       = OPTIONS.wavelet.order;
    J       = OPTIONS.wavelet.nb_levels;
    dt      = diff( OPTIONS.mandatory.DataTime([1 2]) );
    
    [M,N] = size(data);

    % --- definition de l'ondelette
    [wmin,wmax] = Morse_support_spectral(n,m);
    [tmin,tmax] = Morse_support_temporel(n,m);

    % frequence centrale de l'ondelette de Cauchy
    f0          = (((n+0.5)/m).^(1/m))/2/pi; 

    % --- echelles min et max
    a_max = (N-3)*dt/(tmax-tmin);
    a_min = wmax/pi*dt;

    % --- resolution en octaves et en voies
    dj = log2(a_max/a_min)/J;
    aj = a_min*2.^((0:J)*dj);

    % --- coefficients de fourier
    Cwave = repmat(complex(0,0), M,N,J+1);
    WData = [];
    
    % --- options
    OPTIONS.wavelet.scale_analyzed  = fliplr(aj);
    OPTIONS.wavelet.freqs_analyzed  = f0./fliplr(aj);
    OPTIONS.wavelet.Scone  = (f0/dt)./[1E-9,1:((N+1)/2-1),fliplr((1:(N/2-1))),1E-9];
    
    if options_only
        return;
    end
    
 
    % Pad signal - avoid border effects
    sP2 = floor(log2(N));
    pdS = 0;
    if diff([2^sP2 N]) && sP2<11
        pdS = floor( 2^(sP2)-N/2 );
        data = padarray(data, [0 pdS], 'symmetric');
    end


    % --- fft du signal
    f = fft(data,[],2);
    % --- espace de Fourier
    N2 = size(f,2);
    wk = 2.*pi/N2/dt*[1:fix(N2/2)];
    wk = [0., wk, -wk(fix((N2-1)/2):-1:1)];

    % Transform 
    sqrt_aj = sqrt(aj);
    for ai = 1:J+1

        psi_wk = ones(M,1) * ( Morse_w(aj(ai)*wk,n,m) * sqrt_aj(ai) ) ;
        tmpDT = ifft(f.*psi_wk, [], 2);
        Cwave(:,:,J-ai+2) = tmpDT( :, (pdS+1):(end-pdS) );  

    end

    WData = squeeze(Cwave);
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