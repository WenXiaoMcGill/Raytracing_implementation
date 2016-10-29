function [source,newdirec,dist,P,t] = reflection(source,newdirec,plain_coef,L,W,H,P,delta,alpha,c)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% First Wall (North)
temp = dot(newdirec,plain_coef(1,1:3)); % to see if the ray is shot towards the wall
if temp > 0
    % calculate the collision point between the light and the wall
    numer = dot(source,plain_coef(1,1:3))+plain_coef(1,4);
    domi = dot(newdirec, plain_coef(1,1:3));
    point = source - numer/domi * newdirec; 
    % see if the collision point is in the polygon
    if point(1) <= L && point(1) >= 0 && point(3) <= H && point(3) >= 0
        dist = pdist([source; point],'euclidean'); % calculate the distance 
        % calculate the ramaining power of the ray 
        P = P * (1-delta)* exp(-alpha*dist);
        t = dist/c; % calculate the time spent on the ray
        source = point; % set the collision point as the new resource
        % calculate the new direction of the ray
        temp = dot(newdirec,plain_coef(1,1:3))/sqrt(plain_coef(1,1)^2+plain_coef(1,2)^2+plain_coef(1,3)^2);
        newdirec = newdirec-2*temp*plain_coef(1,1:3);
        pos = '1'
        return
    end
end
%% Second wall (South)
temp = dot(newdirec,plain_coef(2,1:3)); % to see if the ray is shot towards the wall
if temp > 0
    % calculate the collision point between the light and the wall
    numer = dot(source,plain_coef(2,1:3))+plain_coef(2,4);
    domi = dot(newdirec, plain_coef(2,1:3));
    point = source - numer/domi * newdirec; 
    % see if the collision point is in the polygon
    if point(1) <= L && point(1) >= 0 && point(3) <= H && point(3) >= 0
        dist = pdist([source; point],'euclidean'); % calculate the distance 
        % calculate the ramaining power of the ray 
        P = P * (1-delta)* exp(-alpha*dist);
        t = dist/c; % calculate the time spent on the ray
        source = point; % set the collision point as the new resource
        % calculate the new direction of the ray
        temp = dot(newdirec,plain_coef(2,1:3))/sqrt(plain_coef(2,1)^2+plain_coef(2,2)^2+plain_coef(2,3)^2);
        newdirec = newdirec-2*temp*plain_coef(2,1:3);
        pos = '2'
        return
    end
end
%% Third wall (East)
temp = dot(newdirec,plain_coef(3,1:3)); % to see if the ray is shot towards the wall
if temp > 0
    % calculate the collision point between the light and the wall
    numer = dot(source,plain_coef(3,1:3))+plain_coef(3,4);
    domi = dot(newdirec, plain_coef(3,1:3));
    point = source - numer/domi * newdirec; 
    % see if the collision point is in the polygon
    if point(2) <= W && point(2) >= 0 && point(3) <= H && point(3) >= 0
        dist = pdist([source; point],'euclidean'); % calculate the distance 
        % calculate the ramaining power of the ray 
        P = P * (1-delta)* exp(-alpha*dist);
        t = dist/c; % calculate the time spent on the ray
        source = point; % set the collision point as the new resource
        % calculate the new direction of the ray
        temp = dot(newdirec,plain_coef(3,1:3))/sqrt(plain_coef(3,1)^2+plain_coef(3,2)^2+plain_coef(3,3)^2);
        newdirec = newdirec-2*temp*plain_coef(3,1:3);
        pos = '3'
        return
    end
end
%% Fourth wall (West)
temp = dot(newdirec,plain_coef(4,1:3)); % to see if the ray is shot towards the wall
if temp > 0
    % calculate the collision point between the light and the wall
    numer = dot(source,plain_coef(4,1:3))+plain_coef(4,4);
    domi = dot(newdirec, plain_coef(4,1:3));
    point = source - numer/domi * newdirec; 
    % see if the collision point is in the polygon
    if point(2) <= W && point(2) >= 0 && point(3) <= H && point(3) >= 0
        dist = pdist([source; point],'euclidean'); % calculate the distance 
        % calculate the ramaining power of the ray 
        P = P * (1-delta)* exp(-alpha*dist);
        t = dist/c; % calculate the time spent on the ray
        source = point; % set the collision point as the new resource
        % calculate the new direction of the ray
        temp = dot(newdirec,plain_coef(4,1:3))/sqrt(plain_coef(4,1)^2+plain_coef(4,2)^2+plain_coef(4,3)^2);
        newdirec = newdirec-2*temp*plain_coef(4,1:3);
        pos = '4'
        return
    end
end
%% Fifth wall (Ceiling)
temp = dot(newdirec,plain_coef(5,1:3)); % to see if the ray is shot towards the wall
if temp > 0
    % calculate the collision point between the light and the wall
    numer = dot(source,plain_coef(5,1:3))+plain_coef(5,4);
    domi = dot(newdirec, plain_coef(5,1:3));
    point = source - numer/domi * newdirec; 
    % see if the collision point is in the polygon
    if point(1) <= L && point(1) >= 0 && point(2) <= W && point(2) >= 0
        dist = pdist([source; point],'euclidean'); % calculate the distance 
        % calculate the ramaining power of the ray 
        P = P * (1-delta)* exp(-alpha*dist);
        t = dist/c; % calculate the time spent on the ray
        source = point; % set the collision point as the new resource
        % calculate the new direction of the ray
        temp = dot(newdirec,plain_coef(5,1:3))/sqrt(plain_coef(5,1)^2+plain_coef(5,2)^2+plain_coef(5,3)^2);
        newdirec = newdirec-2*temp*plain_coef(5,1:3);
        pos = '5'
        return
    end
end
%% Sixth Wall (Floor)
temp = dot(newdirec,plain_coef(6,1:3)); % to see if the ray is shot towards the wall
if temp > 0
    % calculate the collision point between the light and the wall
    numer = dot(source,plain_coef(6,1:3))+plain_coef(6,4);
    domi = dot(newdirec, plain_coef(6,1:3));
    point = source - numer/domi * newdirec; 
    % see if the collision point is in the polygon
    if point(1) <= L && point(1) >= 0 && point(2) <= W && point(2) >= 0
        dist = pdist([source; point],'euclidean'); % calculate the distance 
        % calculate the ramaining power of the ray 
        P = P * (1-delta)* exp(-alpha*dist);
        t = dist/c; % calculate the time spent on the ray
        source = point; % set the collision point as the new resource
        % calculate the new direction of the ray
        temp = dot(newdirec,plain_coef(6,1:3))/sqrt(plain_coef(6,1)^2+plain_coef(6,2)^2+plain_coef(6,3)^2);
        newdirec = newdirec-2*temp*plain_coef(6,1:3);
        pos = '6'
        return
    end
end
end

