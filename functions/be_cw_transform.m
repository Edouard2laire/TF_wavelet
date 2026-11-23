function Cwave = be_cw_transform(data, Psi)
% be_cw_transform: Complex wavelet transform of data (nChan x nTime), using
% the Morse filter bank Psi precomputed by be_cw_init.
%
%   INPUTS:
%       - data : data matrix, size [M x N]
%       - Psi  : filter bank from be_cw_init, size [1 x N x (J+1)]
%
%   OUTPUTS:
%       - Cwave : complex wavelet coefficients, size [M x N x (J+1)]
%
% (Copyright / license header unchanged)

    if size(data,2) ~= size(Psi,2)
        error('be_cw_transform: data has %d time samples but Psi was built for %d. Call be_cw_init again.', ...
              size(data,2), size(Psi,2));
    end

    % Apply all scales at once: M x N x (J+1)
    Cwave = ifft(fft(data, [], 2) .* Psi, [], 2);
end