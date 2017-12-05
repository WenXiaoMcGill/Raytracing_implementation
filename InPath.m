function [flag] = InPath(point,imagepoint,rcv)
% Test if the point is in the path of the two points, if the point is in
% the path, flag =1, else flag = 0

flag = 1;
temp1 = point-imagepoint;
temp2 = point-rcv;
if (temp1(1)*temp2(1)>0||temp1(2)*temp2(2)>0||temp1(3)*temp2(3)>0) || sum(abs(point-imagepoint)<10e-7)==3 || sum(abs(point-rcv)<10e-7)==3
    % This is different from ISM
    flag = 0;
end