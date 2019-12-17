% correlate with base signal
base = alt_filtering(dm34, fs);
%base = alt_filtering(dm54, fs);
TR = {};
%%{
TR{1} = alt_filtering(dm12, fs);
TR{2} = alt_filtering(dm13, fs);
TR{3} = alt_filtering(dm23, fs);
TR{4} = alt_filtering(dm33, fs);
TR{5} = alt_filtering(dm42, fs);
TR{6} = alt_filtering(dm43, fs);
TR{7} = alt_filtering(dm44, fs);
TR{8} = alt_filtering(dm52, fs);
TR{9} = alt_filtering(dm53, fs);
TR{10} = alt_filtering(dm54, fs); % switched with base
%%{
TR{11} = alt_filtering(dm14, fs);
TR{12} = alt_filtering(dm25, fs);
TR{13} = alt_filtering(dm32, fs);
TR{14} = alt_filtering(dm33, fs);
TR{15} = alt_filtering(dm45, fs);
TR{16} = alt_filtering(dm55, fs);
TR{17} = alt_filtering(dm62, fs);
TR{18} = alt_filtering(dm64, fs);
%%{
TR{19} = alt_filtering(dm11, fs);
TR{20} = alt_filtering(dm15, fs);
TR{21} = alt_filtering(dm21, fs);
TR{22} = alt_filtering(dm22, fs);
TR{23} = alt_filtering(dm24, fs);
TR{24} = alt_filtering(dm35, fs);
TR{25} = alt_filtering(dm41, fs);
TR{26} = alt_filtering(dm51, fs);
TR{27} = alt_filtering(dm61, fs);
TR{28} = alt_filtering(dm63, fs);
TR{29} = alt_filtering(dm65, fs);
%}
%{
TR{1} = alt_filtering(dm12, fs);
TR{2} = alt_filtering(dm43, fs);
TR{3} = alt_filtering(dm52, fs);
%}
plt = 0;
n_segs = 64;%ceil(length(base)/fs/0.015);
len_B = length(base)/n_segs;
final = base*0;

for seg = 1:n_segs
    
    sigB = base(1+((seg-1)*len_B):seg*len_B);
    
    [cnt_B, sigB] = remove_interrupts(sigB);
    %sigB(sigB>0.0006) = 0.0006; sigB(sigB<-0.0006) = -0.0006;
    
    sig_out = sigB;
    sig_cnt = cnt_B;
    sig_out(sig_cnt == 0) = 0;
    for k = 1:length(TR)
        
        len_K = length(TR{k})/n_segs;
        sigK = TR{k}(1+((seg-1)*len_K):seg*len_K);
        
        [cnt_K, sigK] = remove_interrupts(sigK);
        %sigK(sigK>0.001) = 0.001; sigK(sigK<-0.001) = -0.001;
        %sigK(sigK>0.0006) = 0.0006; sigK(sigK<-0.0006) = -0.0006;
        
        %sigB(cnt_B == 0) = 0;
        %sigK(cnt_K == 0) = 0;
        
        [cBk, lagBk] = xcorr(sigB,sigK); cBk = cBk/max(cBk);
        [mBk, iBk] = max(cBk); tBk = lagBk(iBk);
         
        sigB(cnt_B == 0) = 0;
        sigK(cnt_K == 0) = 0;
        
        if tBk > 0
            % Base lags behind k
            %{
            sigK = sigK(1:end-tBk);
            cnt_K = cnt_K(1:end-tBk);
            
            if length(sigK) > length(sigB) - tBk
                sigK = sigK(1:length(sigB)-tBk);
                cnt_K = cnt_K(1:length(sigB)-tBk);
            elseif length(sigK) < length(sigB) - tBk
                sigK = [sigK; zeros(length(sigB)-length(sigK)-tBk, 1)];
                cnt_K = [cnt_K; zeros(length(sigB)-length(cnt_K)-tBk, 1)];
            end
            %}
            if length(sigK) > (length(sigB) - tBk)
                sig_cnt(1+tBk:end) = sig_cnt(1+tBk:end) + cnt_K(1:length(sig_cnt)-tBk);
                sig_out(1+tBk:end) = sig_out(1+tBk:end) + sigK(1:length(sig_cnt)-tBk);
            else
                sig_cnt(1+tBk:tBk+length(cnt_K)) = sig_cnt(1+tBk:tBk+length(cnt_K)) + cnt_K;
                sig_out(1+tBk:tBk+length(cnt_K)) = sig_out(1+tBk:tBk+length(cnt_K)) + sigK;
            end
            
            
            if plt == 1
                h = figure(2);
                ax1 = subplot(3,1,1);
                t = 0:1/fs:length(sigB)/fs-1/fs;
                plot(t, sigB)
                
                ax2 = subplot(3,1,2);
                t = 0:1/fs:length(sigK)/fs-1/fs;
                t = t + tBk/fs;
                plot(t, sigK)
                
                ax3 = subplot(3,1,3);
                t = 0:1/fs:length(sig_out)/fs-1/fs;
                plot(t, sig_out)
                
                linkaxes([ax1, ax2, ax3], 'x')
                
                drawnow;
                disp(clock);
                waitfor(h);
                disp(clock);
            end
            
            % can this be extended so that if sigK is long enough
            % it can overlap with more?
            %sig_cnt(1+tBk:end) = sig_cnt(1+tBk:end) + cnt_K;
            %sig_out(1+tBk:end) = sig_out(1+tBk:end) + sigK;
        else
            % Base leads signal K
            %{
            sigK = sigK(1-tBk:end);
            cnt_K = cnt_K(1-tBk:end);
            
            if length(sigK) > length(sigB) + tBk
                sigK = sigK(length(sigK)-length(sigB)-tBk:end);
                cnt_K = cnt_K(length(cnt_K)-length(sigB)-tBk:end);
            elseif length(sigK) < length(sigB) + tBk
                sigK = [zeros(length(sigB)-length(sigK)+tBk, 1); sigK];
                cnt_K = [zeros(length(sigB)-length(cnt_K)+tBk, 1); cnt_K];
            end
            %}
            
            if (length(sigK)-abs(tBk)) > length(sigB)
                
                sig_cnt(1:end) = sig_cnt(1:end) + cnt_K(1+abs(tBk):abs(tBk)+length(sigB));
                sig_out(1:end) = sig_out(1:end) + sigK(1+abs(tBk):abs(tBk)+length(sigB));
            else
                sig_cnt(1:length(cnt_K)-abs(tBk)) = sig_cnt(1:length(cnt_K)-abs(tBk)) + cnt_K(1+abs(tBk):end);
                sig_out(1:length(sigK)-abs(tBk)) = sig_out(1:length(sigK)-abs(tBk)) + sigK(1+abs(tBk):end);
            end
            
            if plt == 1
                h = figure(2);
                ax1 = subplot(3,1,1);
                t = 0:1/fs:length(sigB)/fs-1/fs;
                plot(t, sigB)
                
                ax2 = subplot(3,1,2);
                t = 0:1/fs:length(sigK)/fs-1/fs;
                t = t + tBk/fs;
                plot(t, sigK)
                
                ax3 = subplot(3,1,3);
                t = 0:1/fs:length(sig_out)/fs-1/fs;
                plot(t, sig_out)
                
                linkaxes([ax1, ax2, ax3], 'x')
                
                drawnow;
                disp(clock);
                waitfor(h);
                disp(clock);
            end
            
            %sig_cnt(1:length(sigK)) = sig_cnt(1:length(sigK)) + cnt_K;
            %sig_out(1:length(sigK)) = sig_out(1:length(sigK)) + sigK;
        end
    end
    
    sig_out(sig_cnt == 0) = 0; sig_cnt(sig_cnt == 0) = 1;
    final(1+((seg-1)*len_B):seg*len_B) = sig_out ./ sig_cnt;
        
end
t = 0:1/fs:length(final)/fs-1/fs;
plot(t, final)




