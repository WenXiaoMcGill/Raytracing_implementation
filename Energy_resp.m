function [E] =  Energy_resp(IR, dt, fs)
% This function computes the evergy response of a given impulse response.
% IR: Impulse response vector
% dt: time deviation of each interval
% fs: The sample frequency of the IR

IR_decibel = IR';
% IR_decibel = log10(IR_decibel/(1e-12));
Sample = size(IR_decibel,2);
N = floor(Sample/(dt*fs)); % compute the total number of intervals
E = zeros(1,N);
for n = 1:1:N
    E(n) = sum(IR_decibel((1+floor((n-1)*dt*fs)):Sample));
end
