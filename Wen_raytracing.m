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
L = 50;         % Set up the length of the room
W = 30;         % Set up the width of the room
H = 40;         % Set up the height of the room
S = [10 20 30];    % Set up the position of the sound source
R = [40 10 10];    % Set up the position of the listener
N = 100;         % Set up the number of the rays ommited
%delta = 0.5;    % Set up the absoption coefficient of the walls
%alpha = 0.5;    % Set up the air absorption coefficient
c = 331.4 + 0.6*20; % Set up the sound speed
source = S;
P = 1;   % Set up the power of the sound

% Here I initialize a matrix storing the direction of each ray
ray_direc = zeros(100,3);
ray_direc = ray_direction(ray_direc,N); % Compute all the directions of rays

% input the statistics of the plain
plain_vect = [0 1 0; 0 -1 0; 1 0 0; -1 0 0; 0 0 1; 0 0 -1];
plain_point = [L W H; 0 0 0; L W H; 0 0 0; L W H; 0 0 0];
% calculate the coefficients of each plain Ax+By+Cz+D = 0
row = size(plain_vect,1); % calculate the row of the storing matrix
plain_coef = zeros(row,4);
for n = 1:1:row
    plain_coef(n,1:3) = plain_vect(n,:);
    plain_coef(n,4) = -dot(plain_vect(n,:),plain_point(n,:));
end

%% Reflection computation

% see if the ray is reaching the listener
% calculate the power when the ray reach the listener
% calculate the time when the ray reach the listenner
% restore the power and the time

%% Reflection computation
n = 1;
temp = dot(ray_direc(n,:),plain_coef(1,1:3)); % to see if the ray is shooting towards the wall
if temp > 0
    % calculate the collision point between the light and the wall
    numer = dot(source,plain_coef(1,1:3))+plain_coef(1,4);
    domi = dot(ray_direc(n,:), plain_coef(1,1:3));
    point = source - numer/domi * ray_direc(n,:); 
    % see if the collision point is in the polygon
    if point(1) <= L && point(1) >= 0 && point(3) <= H && point(3) >= 0
        dist = pdist([source; point],'euclidean'); % calculate the distance 
        source = point; % set the collision point as the new resource
        % calculate the new direction of the ray
        temp = dot(ray_direc(n,:)*plain_coef(1,1:3))/sqrt(plain_coef(1,1)^2+plain_coef(1,2)^2+plain_coef(1,3)^2);
        newdirec = ray_direc(n,:)-2*temp*plain_coef(1,1:3);
        % calculate the ramaining power of the ray 
        P = P * exp(alpha*dist) * (1-delta);
        t = dist/c; % calculate the time spent on the ray (remain refined)
    end

end

