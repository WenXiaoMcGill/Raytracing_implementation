function [source,newdirec,dist,P,t] = reflection(source,newdirec,plane,wall,vertex,P,delta,alpha,c)
%  This Function is to find the wall that the ray collides with and compute
%  the new energy and reflecting direction for that ray. I first have to
%  find which wall the ray is colliding with. This is determined by seeing
%  if the ray is shot towards the wall and the crossing point btween the
%  ray and the wall is in the wall's boundary. Once the exact wall is
%  found, I compute the corresponding outputs.

%  Input: source     --- last source position
%         newdirec   --- direction of the emmiting ray
%         plain --- plain coefficients matrix
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


for n = 1:1:size(wall,1)
    % test if the ray is shot towards the wall
    temp = dot(newdirec,plane(n,1:3));
    if temp > 0
        % calculate the collision point between the light and the wall
        numer = dot(source,plane(n,1:3))+plane(n,4);
        domi = dot(newdirec, plane(n,1:3));
        point = source - numer/domi * newdirec;
        % see if the collision point is in the polygon
        % One method used for this study is to form vectors from the point
        % of intersection to each of the vertices on the boundary of the
        % mirror plane. The  cross-products of successive pairs of these
        % vectors are always the vectors that are orthogonal to the
        % reflecting wall. If each of these normal vectors points in the
        % same direction, then the point of intersection is inside the
        % boundary of the reflecting wall. Otherwise it is outside, and
        % such a ray path is invalid.
        comp = [];
        flag = 0;
        for m = 1:1:size(wall,2)
            % ensure the vertices number is not zero
            if wall(n,m) ~= 0 && wall(n,m+1) ~= 0;
                % test if it is the last point to do the last roduct.
                temp = cross((vertex(wall(n,m),:)-point),(vertex(wall(n,m+1),:)-point));
            else if (m == size(wall,2) && wall(n,m)~=0)||(wall(n,m) ~= 0 && wall(n,m+1)==0)
                    temp = cross((vertex(wall(n,m),:)-point),(vertex(wall(n,1),:)-point));
                end
            end
            % check if comparing vector is empty, if it is not empty, do
            % the judgement
            if ~isempty(comp)
                % decide if two vectors are in the same direction
                if temp(1)*comp(1)<0 || temp(2)*comp(2)<0 || temp(3)*comp(3)<0
                    flag = 1;
                    break
                end
            end
            if isempty(comp)
                comp = temp;
            end
        end
        % computation after confirming the the right wall
        if flag ~= 1;
            dist = pdist([source; point],'euclidean'); % calculate the distance
            % calculate the ramaining power of the ray
            P = P * (1-delta)* exp(-alpha*dist);
            t = dist/c; % calculate the time spent on the ray
            source = point; % set the collision point as the new resource
            % calculate the new direction of the ray
            temp = dot(newdirec,plane(n,1:3))/sqrt(plane(n,1)^2+plane(n,2)^2+plane(n,3)^2);
            newdirec = newdirec-2*temp*plane(n,1:3);
            %pos = '1'
            return
        end
    end
    
    %     if point(1) <= L && point(1) >= -1e-14 && point(3) <= H && point(3) >= -1e-14
    %         dist = pdist([source; point],'euclidean'); % calculate the distance
    %         % calculate the ramaining power of the ray
    %         P = P * (1-delta)* exp(-alpha*dist);
    %         t = dist/c; % calculate the time spent on the ray
    %         source = point; % set the collision point as the new resource
    %         % calculate the new direction of the ray
    %         temp = dot(newdirec,plain(1,1:3))/sqrt(plain(1,1)^2+plain(1,2)^2+plain(1,3)^2);
    %         newdirec = newdirec-2*temp*plain(1,1:3);
    %         %pos = '1'
    %         return
    %     end
end
end

