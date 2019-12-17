filename = "trace_KEY_TEST.wav";
trlen = 0.389286;
tr_start = 34.73507;
tr_end = 36.6815;


%trlen = 0.28366;
%filename = "trace_3.wav";
%tr_start = 12.6702;
%tr_end = 14.0895;
%filename = "trace_4.wav";
%tr_start = 7.1988;
%tr_end = 8.6147;
%filename = "trace_5.wav";
%tr_start = 14.6864;
%tr_end = 16.1027;

bp_freq = [1.761*10^6 1.791*10^6];
%bp_freq = [1.175*10^6 1.205*10^6];

% rest should be automated
[raw_full,fs] = audioread(filename);
raw_full = raw_full(tr_start*fs:tr_end*fs);

raw_tr2 = raw_full(1*trlen*fs:2*trlen*fs);
raw_tr3 = raw_full(2*trlen*fs:3*trlen*fs);
raw_tr4 = raw_full(3*trlen*fs:4*trlen*fs);

bp2 = bandpass(raw_tr2, bp_freq, fs);
bp3 = bandpass(raw_tr3, bp_freq, fs);
bp4 = bandpass(raw_tr4, bp_freq, fs);

w = unwrap(angle(hilbert(bp2))); dm2 = diff(w)/(2*pi);
w = unwrap(angle(hilbert(bp3))); dm3 = diff(w)/(2*pi);
w = unwrap(angle(hilbert(bp4))); dm4 = diff(w)/(2*pi);

time_dm = 0:1/fs:length(dm2)/fs-1/fs;

%figure(1)
%subplot(3, 1, 1)
%plot(time_dm, dm2)
%subplot(3, 1, 2)
%plot(time_dm, dm3)
%subplot(3, 1, 3)
%plot(time_dm, dm4)