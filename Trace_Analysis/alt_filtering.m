function y = alt_filtering(dm, fs)
% alternate filtering
%[b_lp, a_lp] = butter(3, 1.2*10^4/(fs/2)); % key test
%dm = dm(0.00225*fs:end);
[b_lp, a_lp] = butter(3, 2.5*10^4/(fs/2)); % unaltered

[b_hp, a_hp] = butter(1, 500/(fs/2), 'high'); % key test
%[b_hp, a_hp] = butter(1, 1.1*10^4/(fs/2), 'high'); % unaltered
y = filter(b_lp, a_lp, dm); y = filter(b_hp, a_hp, y);
y(y>0.005) = 0.005; y(y<-0.005) = -0.005;
%y = smoothdata(y, 'gaussian', 100);
%y = smoothdata(y, 'movmedian', 50);
%plot(time_dm(1:length(y)), y)
end