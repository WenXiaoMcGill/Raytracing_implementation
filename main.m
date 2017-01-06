%% A basic ray-tracing acoustic models in 3D dimensions
% 
% Author: Wen Xiao

%% Initialization
%     ___________________________
%    /|                         /|
%   / |                        / |
%  /__________________________/  |H
%  |  |        * source       |  |
%  |  |_______________________|__|
%  |  /                       | W/
%  | /       listener *       | /
%  |/_________________________|/
%               L
%

clc; clear; 
Initial;
tic
% Here I initialize a matrix storing the direction of each ray
ray_direc = ray_direction(N);
% Set up output storing matrix
out = zeros(N,2);

%% Reflection computation
for n = 1:1:N
source = S;
newdirec = ray_direc(n,:);
P = Po;
Dist = 0;
T = 0;
count = 0;
temp = R-source;
d = norm(cross(temp,newdirec))/norm(newdirec);
while( P > P_thresh)
%while( count < 5)
    % see if the ray is reaching the listener
    temp = R-source; 
    d = norm(cross(temp,newdirec))/norm(newdirec);
    if (d <= r && d ~= 0)
        % calculate the power when the ray reaches the listener
        dtemp1 = sqrt(r^2-d^2);
        dtemp2 = sqrt(norm(temp)^2-d^2);
        w = P * exp(-alpha*dtemp2);
        %Ps = W * 2*dtemp1/(4/3*pi*r^3);
        Ps = P;
        % calculate the time when the ray reaches the listenner
        Ts = T+dtemp2/c;
        % restore the power and the time
        out(n,1) = Ps;
        out(n,2) = Ts;
        break;
    end
    [source,newdirec,dist,P,t] = reflection(source,newdirec,plane,wall,vertex,P,delta,alpha,c);
    T = T+t; % sum the total time
    Dist = Dist + dist; % sum the total distance
    count = count + 1;
end
end

toc
%% visualization
% impulse response
out = sortrows(out,2);
stem(out(:,2),out(:,1));
title('Indoor Ray Tracing Impulse Response')
xlabel('Time (s)')
ylabel('Power')
grid
