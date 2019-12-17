% additional filtering and clipping
LPF = dsp.LowpassFilter('SampleRate', fs, 'FilterType', 'IIR', 'PassbandFrequency', 25*10^3, 'StopbandFrequency', 40*10^3);
HPF = dsp.HighpassFilter('SampleRate', fs, 'FilterType', 'IIR', 'PassbandFrequency', 1000, 'StopbandFrequency', 500);

y2 = step(LPF, dm2); y2 = step(HPF, y2); y2(y2>0.02) = 0.02; y2(y2<-0.02) = -0.02;
y3 = step(LPF, dm3); y3 = step(HPF, y3); y3(y3>0.02) = 0.02; y3(y3<-0.02) = -0.02;
y4 = step(LPF, dm4); y4 = step(HPF, y4); y4(y4>0.02) = 0.02; y4(y4<-0.02) = -0.02;

%figure(2)
%subplot(3, 1, 1)
%plot(time_dm, y2)
%subplot(3, 1, 2)
%plot(time_dm, y3)
%subplot(3, 1, 3)
%plot(time_dm, y4)