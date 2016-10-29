function [ mat ] = ray_direction( mat, N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

d = sqrt(N); % calculate the grids of two coordinates
for n = 0:1:d-1
    %r1 = rand();
    r1 = 0.5;
    for m = 0:1:d-1;
        %r2 = rand();
        r2 = 0.5;
        mat((n*10+m+1),1) = 2*sqrt((n+r1)/d-((n+r1)/d)^2)*cos(2*pi*(m+r2)/d);
        mat((n*10+m+1),2) = 2*sqrt((n+r1)/d-((n+r1)/d)^2)*sin(2*pi*(m+r2)/d);
    end
    mat((n*10+1):(n*10+10),3) = 1-2*(n+r1)/d;
end
end

