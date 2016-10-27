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
%Delta = 0.5;    % Set up the absoption coefficient of the walls
%alpha = 0.5;    % Set up the air absorption coefficient
%c = 331.4 + 0.6T

% Here I initialize a matrix storing the directions of each ray
ray_direc = zeros(100,3);
ray_direc = ray_direction(ray_direc,N); % Compute all the directions of rays
%% Reflection computation



