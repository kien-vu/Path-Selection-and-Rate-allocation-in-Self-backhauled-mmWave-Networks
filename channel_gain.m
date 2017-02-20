function [ gain ] = channel_gain( distance, h )
% This func is used to calculate the channel gain, which taking into
% account small fading, normalized pathloss, interference penalty
%Parameter for LOS at 28 Ghz
alpha3 = 61.4;
beta3 = 2;
interference_penaty = 0.9;
% calculate the noise power 
bandwidth = 1e6; % Hz
boltzmann = 1.381e-23; %constant
kelvin= 290;
noise_figures = 10^(15/10); % 15 dB
noise_power = bandwidth*boltzmann*kelvin*noise_figures; % 1.2665e-13
Normalized_Pass_Loss = alpha3 + beta3 * 10 * log10 (distance) ;%*randn;+ sigma3
Normalized_Pass_Loss = 10.^(-Normalized_Pass_Loss/10);
Normalized_Pass_Loss = 10*log(Normalized_Pass_Loss/noise_power/10);  %*randn;+ sigma3
%%
% d = 40 m, SNR = 3.5751e+03 ~ 58.7916 dB
% d = 43 m, SNR = 791.7179 ~ 43.7162 dB
% d = 43 m, SNR = 791.7179 ~ 43.7162 dB
%%
% gain = interference_penaty * Normalized_Pass_Loss * (rand);
% h = randn + 1i*randn;
gain = interference_penaty * Normalized_Pass_Loss * abs(h);

% h = rand + 1i*rand;
%  abs(h).^2
 
end

