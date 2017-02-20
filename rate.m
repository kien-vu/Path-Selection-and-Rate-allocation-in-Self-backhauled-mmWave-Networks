function [ c ] = rate( a, b  )
% This function is going to calculate the rate for learning game
power = 10; % Assume that each subflow can use maximum transmit power, i.e., 10 wat
[distance_cal, h]  = topology();
    c = ratecal( power, channel_gain( distance_cal(a,b), h(a,b) ) );

end

