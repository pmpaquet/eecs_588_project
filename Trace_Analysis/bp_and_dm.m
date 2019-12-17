function dm = bp_and_dm(y, bp_freq, fs)
    bp = bandpass(y, bp_freq, fs);
    w = unwrap(angle(hilbert(bp)));
    dm = diff(w)/(2*pi);
end