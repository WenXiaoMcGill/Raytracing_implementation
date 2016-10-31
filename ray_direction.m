function [ mat ] = ray_direction( mat, N )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

dn = sqrt(N); % calculate the grids of two coordinates
dn
for n = 0:1:dn-1
    %r1 = rand();
    r1 = 0.5;
    for m = 0:1:dn-1;
        %r2 = rand();
        r2 = 0.5;
        mat((n*dn+m+1),1) = 2*sqrt((n+r1)/dn-((n+r1)/dn)^2)*cos(2*pi*(m+r2)/dn);
        mat((n*dn+m+1),2) = 2*sqrt((n+r1)/dn-((n+r1)/dn)^2)*sin(2*pi*(m+r2)/dn);
    end
    mat((n*dn+1):(n*dn+dn),3) = 1-2*(n+r1)/dn;
end
end

