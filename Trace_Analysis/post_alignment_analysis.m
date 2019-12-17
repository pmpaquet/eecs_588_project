 %x = final(0.00225*fs:end);
 x = final(0.00218*fs:end);
 plot(t(1:length(x)),x)
[pks, locs, w, p] = findpeaks(x, fs);
histogram(pks)
pk_tmp = pks(pks>0.00035); loc_tmp = locs(pks>0.00035);
l1 = loc_tmp(2:end) - loc_tmp(1:end-1);
length(l1(l1>(0.25*w_key)))
l1 = loc_tmp; l1(2:end) = l1(2:end) - l1(1:end-1);
length(l1(l1>(0.25*w_key)))
l1 = loc_tmp(2:end) - loc_tmp(1:end-1);
length(l1(l1>(0.35*w_key)))
length(l1(l1>(0.5*w_key)))
length(l1(l1>(0.75*w_key)))
length(l1(l1>(0.85*w_key)))
length(l1(l1>(0.95*w_key)))
length(l1(l1>(1.2*w_key)))
length(l1(l1>(1*w_key)))
length(l1(l1>(0.75*w_key)))
bit1_indexer = [1;l1>0.75*w_key];
sum(bit1_indexer)
 