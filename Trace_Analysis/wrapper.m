% all together now
run('get_signals.m');

y2 = alt_filtering(dm2, fs);
y3 = alt_filtering(dm3, fs);
y4 = alt_filtering(dm4, fs);

run('avg_signal.m');

plot(time_dm(1:length(new_sig)), new_sig);