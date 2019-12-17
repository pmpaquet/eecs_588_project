% get secret key
w_key = 1/(2036*fs/length(x));
prev_loc = 0;
secret_key = "";
for i = 1:length(bit1_locations)
    num_zeros = round((bit1_locations(i) - prev_loc)/(w_key)) - 1;
    for j = 1:num_zeros
        secret_key = secret_key + "0";
    end
    secret_key = secret_key + "1";
    prev_loc = bit1_locations(i);
end

num_zeros = round(((length(x)/fs) - prev_loc)/w_key) - 1;
for j = 1:num_zeros
    secret_key = secret_key + "0";
end
