function [ mat_x mat_y ] = Img_pos( Dim,R,L,W )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

Center_x = (Dim+1)/2; % locate x coordinate of the center case
Center_y = (Dim+1)/2; % locate y coordinate the center case
mat_x(Center_x,Center_y) = R(1); % save the original positon of the listener
mat_y(Center_x,Center_y) = R(2); % save the original positon of the listener

end

