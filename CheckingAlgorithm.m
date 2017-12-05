Debug = 0;
Start_ray = 1;
Check_plane = 1;

%%------------------------------------------------------------------------
Band_num = 4;

% Initial_MMR
Initial_Pollack
% Initial_Tanna
% Initial_Wirth

% Decide the display of the wall filters
DispFilter = 0;
[B,A] = BPFilter(f(Band_num), fs, 2, DispFilter);

% choose file type
load (addr1a) % type the correct address acording to the Workspace
y = fftdeconv(response, x, 0.99); % compute the measured IR

% choose source and reveiver positions
position = 1;
% choose the reference direction
direc = 1; 
src = src(position,:); % Speaker postion
rcv = rcv(position,:); % Mic position
ref_direc = ref_direc(direc,:);

% choose display type
overlap_only = 1; % only display overlaped IR or not 1,0,-1
imprange = [4]; % set the second order impulse that we want to see 0==all -1==none
refl_order = 2; % reflection order 1,2,-1

N = 210*210;            % Set up the number of the rays emitted must be power of 2
% r = log10(V)*norm(src-rcv)*sqrt(4/N);    % Set up the radius of listener
r = 0.4;
delta = 0.005;        % Set up the absoption coefficient of the walls
alpha = 0.001;        % Set up the air absorption coefficient
c = 331.4 + 0.6*23.6; % Set up the sound speed: c + 0.6T, where T is room temperature
Po = 1;               % Set up the power of the sound
P_thresh = 0.001*Po;  % Set up the threshhold of the ray energy
Time = 0.1;               % Time 
Fs = fs;               % Sample Rate

% generate the function of every wall
plane = TPlane(wall,wnum,vertex);

% Here I initialize a matrix storing the direction of each ray
% ray_direc = ray_direction(N);

% Set up output storing matrix
out = zeros(N,2+refl_order);

% Room frame display %----------------------------------------------------%
if Debug == 1
    figure;
    xx = zeros(wnum*size(wall,2),2);
    yy = zeros(wnum*size(wall,2),2);
    zz = zeros(wnum*size(wall,2),2);
    for u = 1:1:wnum % wall number counter
        for v = 1:1:size(wall,2) % vertice counter for every wall
            if (v == size(wall,2) && wall(u,v)~=0)||(wall(u,v) ~= 0 && wall(u,v+1)==0)
                xx((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),1),vertex(wall(u,1),1)];
                yy((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),2),vertex(wall(u,1),2)];
                zz((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),3),vertex(wall(u,1),3)];
                continue
            else if v < size(wall,2) && wall(u,v+1)~=0
                    xx((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),1),vertex(wall(u,v+1),1)];
                    yy((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),2),vertex(wall(u,v+1),2)];
                    zz((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),3),vertex(wall(u,v+1),3)];
                end
            end
        end
    end
    xx = xx';
    yy = yy';
    zz = zz';
    plot3(xx, yy, zz,'k-','LineWidth',2)
    hold on
    % plot the source and the receiver position
    plot3(src(1,1),src(1,2),src(1,3),'.','MarkerSize',50) % Source position
    plot3(rcv(1,1),rcv(1,2),rcv(1,3),'.','MarkerSize',50) % Receiver position
end
%------------------------------------------------------------------------%
%% Reflection computation
for n = Start_ray:1:N
    source = src(1,:);
    tpoint = source;
    newdirec = ray_direc(n,:);
    P = Po;
    Dist = 0;
    T = 0;
    count = 0;
    while(count <= refl_order)
        % check if the ray is reaching the listener
        temp = rcv(1,:)-source;
        d = norm(cross(temp,newdirec))/norm(newdirec);
        % compute the direct sound path
%       d_direc = norm(source-rcv);
%       E_direc = P* exp(-alpha*d_direc)/(4*pi*d_direc^2);
%       T_direc = d_direc/c;
        if (d <= r && d ~= 0)
            % calculate the power when the ray reaches the listener
            dtemp1 = sqrt(r^2-d^2);
            dtemp2 = sqrt(norm(temp)^2-d^2);
            spk_coef = direc_judge(ref_direc,ray_direc(n,:));
            % calculate the time when the ray reaches the listenner
            Ts = T+dtemp2/c;
            % restore the power and the time
%             out(n,1) = P* exp(-alpha*dtemp2)*spk_coef/(4*pi*(c*Ts)^2)*2*dtemp1/(4/3*pi*r^3);
            out(n,1) = P* exp(-alpha*dtemp2)*spk_coef*2*dtemp1/(4/3*pi*r^3);
            out(n,2) = Ts;
            break
        end
        [source,newdirec,dist,P,t,out] = refl(source,newdirec,plane,wall,vertex,P,alpha,c,out,n,count,beta);
%         [source,newdirec,dist,P,t] = refl(source,newdirec,plane,wall,vertex,P,beta,alpha,c,f_concrete,f_wood,f_lime);
        % Ray path display-------------------------------------------------%
        if Debug == 1
            disp([source; newdirec]);
            plot3(source(1),source(2),source(3),'.','MarkerSize',5) % Source position
            plot3([tpoint(1) source(1)],[tpoint(2) source(2)],[tpoint(3) source(3)],'Color',[0 0.003*count 0],'LineStyle','-','LineWidth',2)
            plot3([source(1) source(1)+newdirec(1)],[source(2) source(2)+newdirec(2)],[source(3) source(3)+newdirec(3)],'Color',[0 0.003*count 0],'LineStyle','-','LineWidth',2)
            tpoint = source;
        end
        %------------------------------------------------------------------%
        T = T+t; % sum the total time
        if T > Time
            break
        end
        Dist = Dist + dist; % sum the total distance
        count = count + 1;
    end
end