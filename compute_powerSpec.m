function powerSpec = compute_powerSpec( OPL, k0, rep )
%COMPUTE_POWERSPEC given OPL and k0, calculate power spectrum of a photon packet

        phase = OPL' * k0;                                % phase of an E-field
        phase_rep = repmat(phase,[rep,1]);         % nPeriod = simulated time length * fsampling = 4us*25.25MHz = 101
        E_tag = exp(1i * phase_rep);                  % E-field of light
        FT_Etag = fftshift(fft(E_tag));                  % Fourier Transform of E field
        powerSpec = abs(FT_Etag).^2;                       % Power Spectrum of E field
        powerSpec = powerSpec./sum(powerSpec);                   % Normalized power spectrum

end

