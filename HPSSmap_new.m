function [ harm, perc ] = HPSSmap( spect, iterations )
% HPSS - MAXIMUM A POSTERIORI
% Adam Szaruga, 7/31/2014
%
% Modified some code for accelaration, not changing input and output.
% Daniel Nango, 6/25/2016
%
% Separates a spectrogram into its harmonic and percussive elements
%
% #ARGUMENTS#
% spect - the input spectrogram
% iterations - number of iterations to run
%
% Returns two spectrograms: 
% harm - The spectrogram with supressed percussive elements
% perc - The spectrogram with supressed harmonic elements
%
% #USAGE#
% [s, Fs] = audioread( path );
% spect = Spect(s);
% [harm, perc] = HPSSmap(spect, 10);
% h = iSpect(harm);
% p = iSpect(perc);
% hPlayer = audioplayer(h, Fs);
% pPlayer = audioplayer(p, Fs);
% play(hPlayer);
% play(pPlayer);
%

% Extract magnitude and phase from the spectrogram
mag = abs(spect);
phase = angle(spect);

% initialize harmonic and percussive spectrograms
harm = mag./2;
perc = mag./2;

% check for scenario where harm and perc are zero, equalize
zero_harm = (harm == 0);
zero_perc = (perc == 0);
harm(zero_harm & zero_perc) = 0.5;
perc(zero_harm & zero_perc) = 0.5;

%initialize auxillary variables
mh = harm./(harm + perc);
mp = perc./(harm + perc);

for n = 1:iterations 
    % calculate variances
    varh = var(harm(:));
    varp = var(perc(:));
    
    % calculate prelim variables
    ah = (2/varh) + 2;
    bh = bhCalc(harm, varh); 
    ch = mh.*mag.*2;
    
    ap = (2/varp) + 2;
    bp = bpCalc(perc, varp); 
    cp = mp.*mag.*2;
    
    % update harmonic and percussive spectrograms
    harm = ((bh + sqrt(bh.^2 + ch.*(4*ah)))./(2*ah)).^2;
    perc = ((bp + sqrt(bp.^2 + cp.*(4*ap)))./(2*ap)).^2;
    
    % update auxillary variables
    mh = harm./(harm + perc);
    mp = perc./(harm + perc);
end

% reapply phase to magnitude
harm = harm.*exp(1i*phase);
perc = perc.*exp(1i*phase);

end


function bh = bhCalc(harm, varh)

harm = sqrt(harm);

bh = zeros(size(harm,1), size(harm,2));

bh(:,1) = harm(:,2);

for k = 2:size(bh,2)-1
    bh(:,k) = harm(:, k-1) + harm(:, k+1);
end

bh(:,size(bh,2)) = harm(:,size(bh,2)-1);

bh = bh./varh;

end


function bp = bpCalc(perc, varp)

perc = perc';
perc = sqrt(perc);

bp = zeros(size(perc,1), size(perc,2));


bp(:,1) = perc(:,2);

for n = 2:size(bp,2)-1
    bp(:,n) = perc(:,n-1) + perc(:,n+1);
end

bp(:,size(bp,2)) = perc(:,size(bp,2)-1);
bp = bp';

bp = bp./varp;
end
