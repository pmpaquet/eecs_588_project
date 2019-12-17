% correlate with base signal
%base = alt_filtering(dm34, fs);
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
TR{10} = alt_filtering(dm54, fs);
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
TR{30} = alt_filtering(dm34, fs);
%}
%{
TR{1} = alt_filtering(dm12, fs);
TR{2} = alt_filtering(dm43, fs);
TR{3} = alt_filtering(dm52, fs);
%}
plt = 0;
n_segs = 64;%ceil(length(base)/fs/0.015);

new_segs = {};

for seg = 1:n_segs
    
    %sigB = base(1+((seg-1)*len_B):seg*len_B);
    
    %st_ind = 1;
    
    %[cnt_B, sigB] = remove_interrupts(sigB);
    %sigB(sigB>0.0006) = 0.0006; sigB(sigB<-0.0006) = -0.0006;
    
    %sig_out = sigB;
    %sig_cnt = cnt_B;
    %sig_out(sig_cnt == 0) = 0;
    for k = 2:length(TR)
        
        len_B = length(TR{k-1})/n_segs;
        sigB = TR{k-1}(1+((seg-1)*len_B):seg*len_B);
        
        len_K = length(TR{k})/n_segs;
        sigK = TR{k}(1+((seg-1)*len_K):seg*len_K);
        
        [cnt_B, sigB] = remove_interrupts(sigB);
        [cnt_K, sigK] = remove_interrupts(sigK);
        
        if k == 2
            sig_out = zeros(length(sigB), 1); 
            sig_cnt = zeros(length(cnt_B), 1);
            ind_B = 1;
        end
        
        [cBk, lagBk] = xcorr(sigB,sigK); cBk = cBk/max(cBk);
        [mBk, iBk] = max(cBk); tBk = lagBk(iBk);
         
        sigB(cnt_B == 0) = 0;
        sigK(cnt_K == 0) = 0;
        
        %save where base(1) exists
        base_ind = 1 + max(-tBk, 0);
        
        [tmp_cnt, tmp_out] = align_BK(cnt_B, sigB, cnt_K, sigK, tBk);
        
        % add together
        
        if k == 2
            sig_cnt = tmp_cnt;
            sig_out = tmp_out;
        else
            %sig_out(sig_cnt == 0) = 0; sig_cnt(sig_cnt == 0) = 1;
            [sig_cnt, sig_out] = align_BK(sig_cnt, sig_out, cnt_K, sigK, prev_base_ind - base_ind);
        end
        prev_base_ind = base_ind;
        
        
    end
    
    sig_out(sig_cnt == 0) = 0; sig_cnt(sig_cnt == 0) = 1;
    new_segs{seg} = sig_out ./ sig_cnt;
        
end

for i = 1:n_segs
    % In theory there shouldn't be overlap between 3 segments
   if i == 1
       final = new_segs{i};
   else
       sigB = new_segs{i-1};
       sigK = new_segs{i};
       
       [cBk, lagBk] = xcorr(sigB,sigK); cBk = cBk/max(cBk);
       [mBk, iBk] = max(cBk); tBk = lagBk(iBk);
       
       if tBk < 0
           
           final = [final; sigK];
           continue
       elseif tBk == 0
           final = [final; sigK];
           continue
       end
       
       ind_end = length(cnt_B) + max(length(cnt_K) + tBk - length(cnt_B), 0);
       if ind_end < length(sigB)
           
           break
       end
       
       ind_start = ind_end - length(sigK) + 1;
       offset = length(final) - length(sigB);
       
       final(ind_start+offset:end) = (final(ind_start+offset:end) + sigK(1:length(sigB)-tBk))/2;
       final = [final; sigK(length(sigB)-tBk+1:end)];
   
       %[final_cnt, final] = align_BK(ones(length(final),1), final, ones(length(new_segs{i}),1), new_segs{i}, tFinSeg);
       %final = final ./ final_cnt;
   end
       
end

t = 0:1/fs:length(final)/fs-1/fs;
plot(t,final)
    

t = 0:1/fs:length(final)/fs-1/fs;
plot(t, final)
