function [source,newdirec,dist,P,t,out] = refl(source,newdirec,plane,wall,vertex,P,alpha,c,out,nn,count,beta)
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
    newdirec = newdirec/norm(newdirec);
    temp = dot(newdirec,plane(n,1:3));
    if temp > 0
        % calculate the collision point between the light and the wall
        numer = dot(source,plane(n,1:3))+plane(n,4);
        domi = dot(newdirec, plane(n,1:3));
        point = source - numer/domi * newdirec;
        % debugging display -------------------------------------------%
%         if n == 9
%             disp('The intersection point is');
%             disp([point]);
%         end
        %---------------------------------------------------------------%
        %(1) test if the ray is in the opposite direction
        temp = point-source;
        if dot(temp,newdirec)<0

            continue
        end
        %(2) see if the collision point is in the polygon
        flag1 = 1;
        flag2 = 0;
        flag1 = BoundaryJudge(wall,vertex,point,n);
        % debugging display -------------------------------------------%
%         if n == 9
%             disp('In the polygon?');
%             disp(flag1);
%         end
        %---------------------------------------------------------------%
        if flag1 == 0
            continue
        end
        %(3) test if there exist any obstructions in the path
        for m = 1:1:size(wall,1)
            if m ~= n && dot(point-source,plane(m,1:3))>0
                tpoint = CrossPoint(point,source,plane,m);
                temp = BoundaryJudge(wall,vertex,tpoint,m);
                if temp == 1
                    flag2 = InPath(tpoint,point,source);
                    if flag2 == 1
                        % debugging display -------------------------------------------%
%                         if n == 9
%                             disp('Path obstruction wall is');
%                             disp([m]);
%                         end
                        %--------------------------------------------------------------%
                        break
                    end
                end
            end
        end

        % computation after confirming the the right wall
        if flag1 == 1 && flag2 ~= 1
            dist = pdist([source; point],'euclidean'); % calculate the distance
            % calculate the ramaining power of the ray
            P = P * exp(-alpha*dist)*beta(n);
            out(nn,2+count+1) = n;
            t = dist/c; % calculate the time spent on the ray
            % calculate the new direction of the ray after reflection
            imgsource = Mirror(wall,vertex,source,n,plane);
            newdirec = point - imgsource;
            source = point; % set the collision point as the new resource
            % remove the deviation during the calculation
            if isempty(find(abs(source)<1e-15, 1)) == 0
                source(find(abs(source)<1e-15,1)) = 0;
            end
%             temp = dot(newdirec,plane(n,1:3))/norm(plane(1:3));
%             newdirec = newdirec-2*temp*plane(n,1:3);
            return
        end
    end
end
