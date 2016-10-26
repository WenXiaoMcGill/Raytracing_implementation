%% A basic ray-tracing acoustic models in 2D dimensions
% Here I use the rectangular room. 
% Reference: 
% Author: Wen Xiao

%% Initialization
%  The model 
%   __________________________
%  |                          |
%  |   * source               |
%  |                          | W
%  |         listener *       |
%  |__________________________|
%               L
%
L = 50;         % Set up the length of the room
W = 30;         % Set up the width of the room
S = [10 20];    % Set up the position of the sound source
R = [40 10];    % Set up the position of the listener
%Delta = 0.5;    % Set up the absoption coefficient of the walls
%c = 331.4 + 0.6T

% Here I initialize a matrix storing the image positions of the listener
% The number of the dimention should be odd to better locate the center
% case

Dim = 3;              % the dimention of the image matrix (Dim*Dim). 
Img_lis_x = zeros(Dim,Dim); % Set up the image matrix in x coordinate
Img_lis_y = zeros(Dim,Dim); % Set up the image matrix in y coordinate

% Compute all the positions of image listeners
[Img_lis_x,Img_lis_y] = Img_pos(Dim,R,L,W);

%%



