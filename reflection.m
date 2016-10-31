function [source,newdirec,dist,P,t] = reflection(source,newdirec,plain_coef,L,W,H,P,delta,alpha,c)
%  This Function is to find the wall that the ray collides with and compute
%  the new energy and reflecting direction for that ray. I first have to
%  find which wall the ray is colliding with. This is determined by seeing
%  if the ray is shot towards the wall and the crossing point btween the
%  ray and the wall is in the wall's boundary. Once the exact wall is
%  found, I compute the corresponding outputs.

%  Input: source     --- last source position
%         newdirec   --- direction of the emmiting ray
%         plain_coef --- plain coefficients matrix
%         L          --- the length of the room
%         W          --- the width of the room
%         H          --- the height of the room
%         P          --- the power of the sound
%         delta      --- the absoption coefficient of the walls
%         alpha      --- the air absorption coefficient
%         c          --- velocity of sound
%  Output: source    --- new source position
%          newdirec  --- new direction of the ray after reflection
%          dist      --- distance form the source to the collition point
%          P         --- power remained after collision
%          t         --- time spent for the refelction

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
    point
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
    % dot(point,plain_coef(5,1:3))+plain_coef(5,4) %test
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
        %dot(cross(olddirec,newdirec),plain_coef(5,1:3)) % test
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

