function output = FWHM(waveform, negative)
if negative==1

    minV = min(waveform);
    output = numel(waveform(waveform<=0.5*minV));
end
end