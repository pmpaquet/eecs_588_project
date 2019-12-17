%align signals
[c32, lag32] = xcorr(s3,s2); c32 = c32/max(c32);
[c43, lag43] = xcorr(s4,s3); c43 = c43/max(c43);
[c42, lag42] = xcorr(s4,s2); c42 = c42/max(c42);

[m32, i32] = max(c32); t32 = lag32(i32);
[m43, i43] = max(c43); t43 = lag43(i43);
[m42, i42] = max(c42); t42 = lag42(i42);

c32 = 0; c43 = 0; c42 = 0;
lag32 = 0; lag43 = 0; lag42 = 0;


% 2 [-------+++++---]
% 3 [-----+++++-----]
% 4 [-+++++---------]

% 2 [-----+++++---]
% 3 [-----+++++-----]
% 4     [-+++++-----]

if t32 > 0 && t42 > 0
    first = s2;
    tST = -abs(t43);
    if t43 > 0
        second = s3; third = s4;
        tFS = -abs(t32); tFT = -abs(t42);
    else
        second = s4; third = s3;
        tFS = -abs(t42); tFT = -abs(t32);
    end
elseif t42 < 0 && t43 < 0
    first = s4;
    tST = -abs(t32);
    if t32 > 0
        second = s2; third = s3;
        tFS = -abs(t42); tFT = -abs(t43);
    else
        second = s3; third = s2;
        tFS = -abs(t43); tFT = -abs(t42);
    end
else
    first = s3;
    tST = -abs(t42);
    if t42 > 0
        second = s2; third = s4;
        tFS = -abs(t32); tFT = -abs(t43);
    else
        second = s4; third = s2;
        tFS = -abs(t43); tFT = -abs(t32);
    end
end

% Remove interupts
cnt_first = remove_interrupts(first);
cnt_second = remove_interrupts(second);
cnt_third = remove_interrupts(third);

sig_cnt = cnt_second;
sig_out = second;
sig_out(sig_cnt == 0) = 0;

% first
first = first(1:end+tFS);
cnt_first = cnt_first(1:end+tFS);
tmp = first; tmp(cnt_first == 0) = 0;
sig_cnt(1-tFS:end) = sig_cnt(1-tFS:end) + cnt_first;
sig_out(1-tFS:end) = sig_out(1-tFS:end) + tmp;

% third
third = third(1-tST:end);
cnt_third = cnt_third(1-tST:end);
tmp = third; tmp(cnt_third == 0) = 0;
sig_cnt(1:length(third)) = sig_cnt(1:length(third)) + cnt_third;
sig_out(1:length(third)) = sig_out(1:length(third)) + tmp;

sig_out(sig_cnt == 0) = 0; sig_cnt(sig_cnt == 0) = 1;

% output
sig_out = sig_out ./ sig_cnt;

%{
h = figure(3);
subplot(4,1,1)
y = zeros(length(second), 1); y(1-tFS:end) = first;
plot(y)
plot(time_dm(1:length(y)), y)
subplot(4, 1, 2)
plot(time_dm(1:length(y)), second)
subplot(4,1,3)
y = zeros(length(second), 1); y(1:length(third)) = third;
plot(time_dm(1:length(y)), y)
subplot(4, 1, 4)
plot(time_dm(1:length(y)), sig_out)

drawnow;
disp(clock);
waitfor(h);
disp(clock);
%}



