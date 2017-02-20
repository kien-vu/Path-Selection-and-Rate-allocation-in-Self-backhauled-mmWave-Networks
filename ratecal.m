function [ R ] = ratecal( power, gain )
%   Detailed explanation goes here
%     gain = channel_gain(80,0.1);
%     power = 5;
%     type = 0;
        R = log(1 + power * gain.^2);

end

