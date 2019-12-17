function [cnt, arr] = remove_interrupts(arr)
    cnt = ones(length(arr), 1);
    spread=1000;
    rand_spread = 0.00075;
    for i = 1:length(arr)
        if abs(arr(i)) >= 0.001%25
            if i > spread && (length(arr) - i) > spread
                cnt(i-spread:i+spread) = 0;
                arr(i-spread:i+spread) = (rand(2*spread+1, 1)-0.5)*rand_spread*2;
            elseif i <= spread
                cnt(1:i+spread) = 0;
                arr(1:i+spread) = (rand(i+spread, 1)-0.5)*rand_spread*2;
            else
                cnt(i-spread:end) = 0;
                arr(i-spread:end) = (rand(length(arr(i-spread:end)), 1)-0.5)*rand_spread*2;
            end
        end
    end
    %arr(cnt == 0) = 0;                
end