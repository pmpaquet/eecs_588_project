n_segs = 6;
seg_ln = length(y2)/n_segs;

new_sig = zeros(length(y2),1);
for seg_num = 1:n_segs
s2 = y2(1+(seg_num-1)*seg_ln:seg_num*seg_ln);
s3 = y3(1+(seg_num-1)*seg_ln:seg_num*seg_ln);
s4 = y4(1+(seg_num-1)*seg_ln:seg_num*seg_ln);
%%{
s2(s2>0.001) = 0.001; s2(s2<-0.001) = -0.001;
s3(s3>0.001) = 0.001; s3(s3<-0.001) = -0.001;
s4(s4>0.001) = 0.001; s4(s4<-0.001) = -0.001;
%}

%{
s2(abs(s2)>0.001) = 0;
s3(abs(s3)>0.001) = 0;
s4(abs(s4)>0.001) = 0;
%}
run('align_tr.m');
new_sig(1+(seg_num-1)*seg_ln:seg_num*seg_ln) = sig_out;
end