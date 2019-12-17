function [tmp_cnt, tmp_out] = align_BK(cnt_B, sigB, cnt_K, sigK, tBk)
    
    if tBk > 0
            % Base lags behind k
            % Base index is 1
            ind_end = length(cnt_B) + max(length(cnt_K) + tBk - length(cnt_B), 0);
            
            tmp_cnt = zeros(ind_end, 1);
            tmp_out = zeros(ind_end, 1);
                
            tmp_cnt(1:length(cnt_B)) = cnt_B;
            tmp_out(1:length(sigB)) = sigB;
            
            %[length(tmp_cnt(1+tBk:end)), length(cnt_K), tBk];
                
            tmp_cnt(1+tBk:tBk+length(cnt_K)) = tmp_cnt(1+tBk:tBk+length(cnt_K)) + cnt_K;
            tmp_out(1+tBk:tBk+length(sigK)) = tmp_out(1+tBk:tBk+length(sigK)) + sigK;
    else
            % Base leads signal K
            % base index is abs(tBk)
            ind_end = length(sigB) + max(length(sigK)-abs(tBk)-length(sigB), 0);
            
            tmp_cnt = zeros(ind_end + abs(tBk)+1, 1);
            tmp_out = zeros(ind_end + abs(tBk)+1, 1);
                
            tmp_cnt(1+abs(tBk):abs(tBk)+length(cnt_B)) = cnt_B;
            tmp_out(1+abs(tBk):abs(tBk)+length(sigB)) = sigB;
                
            tmp_cnt(1:length(cnt_K)) = tmp_cnt(1:length(cnt_K)) + cnt_K;
            tmp_out(1:length(cnt_K)) = tmp_out(1:length(cnt_K)) + sigK;
    end
end