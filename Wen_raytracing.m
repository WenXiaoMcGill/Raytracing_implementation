%% A basic ray-tracing acoustic models in 3D dimensions
% Here I use the rectangular room. 
% Reference: 
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

L = 50;             % Set up the length of the room
W = 30;             % Set up the width of the room
H = 40;             % Set up the height of the room
S = [10 20 30];     % Set up the position of the sound source
R = [40 10 10];     % Set up the position of the listener
r = 2.5;             % Set up the radius of listener
N = 20*20;            % Set up the number of the rays ommited must be power of 2
delta = 0.005;        % Set up the absoption coefficient of the walls
alpha = 0.001;        % Set up the air absorption coefficient
c = 331.4 + 0.6*20; % Set up the sound speed
Po = 1/N;            % Set up the power of the sound
P_thresh = 0.0001*Po;   % Set up the threshhold of the ray energy

% Here I initialize a matrix storing the direction of each ray
ray_direc = zeros(N,3);
ray_direc = ray_direction(ray_direc,N); % Compute all the directions of rays

% input the statistics of the plain
plain_vect = [0 1 0; 0 -1 0; 1 0 0; -1 0 0; 0 0 1; 0 0 -1];
plain_point = [L W H; 0 0 0; L W H; 0 0 0; L W H; 0 0 0];
% calculate the coefficients of each plain Ax+By+Cz+D = 0
row = size(plain_vect,1); % calculate the row of the storing matrix
plain_coef = zeros(row,4);
for n = 1:1:row
    plain_coef(n,1:3) = plain_vect(n,:); % denote A B C
    plain_coef(n,4) = -dot(plain_vect(n,:),plain_point(n,:)); % calculate D
end
out = zeros(N,2); % Set up output storing matrix
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
    % see if the ray is reaching the listener
    temp = R-source; 
    d = norm(cross(temp,newdirec))/norm(newdirec);
    if (d <= r)
        % calculate the power when the ray reach the listener
        dtemp1 = sqrt(r^2-d^2);
        dtemp2 = sqrt(norm(temp)^2-d^2);
        w = P * exp(-alpha*dtemp2);
        %Ps = W * 2*dtemp1/(4/3*pi*r^3);
        Ps = P;
        % calculate the time when the ray reach the listenner
        Ts = T+dtemp2/c;
        % restore the power and the time
        out(n,1) = Ps;
        out(n,2) = Ts;
        break;
    end
    [source,newdirec,dist,P,t] = reflection(source,newdirec,plain_coef,L,W,H,P,delta,alpha,c);
    T = T+t; % sum the total time
    Dist = Dist + dist; % sum the total distance
    count = count + 1;
end
end
%% visualization
% impulse response
out = sortrows(out,2);
stem(out(:,2),out(:,1));


