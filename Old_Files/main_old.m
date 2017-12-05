%% A basic ray-tracing acoustic models in 3D dimensions
% 
% Author: Wen Xiao

%% Initialization
%     ___________________________
%    /|                         /|
%   / |                        / |
%  /__________________________/  |H
%  |  |        * source       |  |
%  |  |_______________________|__|
%  |  /                       | W/
%  | /       listener *       | /
%  |/_________________________|/
%               L

clc; clear;
tic % Begin calculating the computation time

% Initialize the geometry data of rooms
% Initial_Rect
% Initial_Queen_Marry_classroom
Initial_MMR

position = 1; % choose source and reveiver positions
diplay_audio = 'mmr1.wav';
overlap_only = 1; % only display overlaped IR or not 1,0,-1
imprange = [4]; % set the second order impulse that we want to see 0==all -1==none
refl_order = 2; % reflection order 1,2,-1
wallrange = [1:wnum]; % walls that waves will hit in the second reflection

r = 0.4;              % Set up the radius of listener
N = 230*230;            % Set up the number of the rays ommited must be power of 2
delta = 0.005;        % Set up the absoption coefficient of the walls
alpha = 0.001;        % Set up the air absorption coefficient
c = 331.4 + 0.6*23.6; % Set up the sound speed: c + 0.6T, where T is room temperature
Po = 1/1500;            % Set up the power of the sound
P_thresh = 0.001*Po;   % Set up the threshhold of the ray energy
Time = 0.1;               % Time 
Fs = 48000;            % Sample Rate

% Set up the position of the sound source and the listener
if position == 1
    % mmr1
    src = [7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser]; % MMR Speaker
    rcv = [8.42+dLaser, 15.14+dLaser, 3.55+0.13+dLaser]; % MMR Mic
else if position == 2
        % mmr2
        src = [7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser]; % MMR Speaker
        rcv = [11+dLaser, 17.67-dLaser, 3.55+0.13+dLaser]; % MMR Mic
    else if position == 3
            % mmr3
            src = [7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser]; % MMR Speaker
            rcv = [13.6780-dLaser, 15.86-dLaser, 3.55+0.13+dLaser]; % MMR Mic
        else if position == 4
                % mmr4
                src = [7.53+dLaser, 7.94+dLaser, 3.63+0.15+dLaser]; % MMR Speaker
                rcv = [8.71+dLaser, 16.76-dLaser, 3.55+0.13+dLaser]; % MMR Mic
            end
        end
    end
end

% generate the function of every wall
plane = TPlane(wall,wnum,vertex);

% Here I initialize a matrix storing the direction of each ray
ray_direc = ray_direction(N);
% Set up output storing matrix
out = zeros(N,2);

%% Reflection computation
for n = 1:1:N
    source = src;
    newdirec = ray_direc(n,:);
    P = Po;
    Dist = 0;
    T = 0;
    count = 0;
    temp = rcv-source;
    d = norm(cross(temp,newdirec))/norm(newdirec);
    tpoint = src;
%     while( P > P_thresh && count <= 2)
    while( count <= 3)
        % see if the ray is reaching the listener
        temp = rcv-source;
        d = norm(cross(temp,newdirec))/norm(newdirec);
        if (d <= r && d ~= 0)
            % calculate the power when the ray reaches the listener
            dtemp1 = sqrt(r^2-d^2);
            dtemp2 = sqrt(norm(temp)^2-d^2);
            spk_coef = direc_judge(src,rcv,ray_direc(n,:));
            w = P * exp(-alpha*dtemp2);
            %Ps = W * 2*dtemp1/(4/3*pi*r^3);
            Ps = P;
            % calculate the time when the ray reaches the listenner
            Ts = T+dtemp2/c;
            % restore the power and the time
            out(n,1) = Ps;
            out(n,2) = Ts;
            break;
        end
        [source,newdirec,dist,P,t] = refl(source,newdirec,plane,wall,vertex,P,beta,alpha,c);
        %------------------------------------------------------------------%
%         plot3(source(1),source(2),source(3),'.','MarkerSize',20) % Source position
%         plot3([tpoint(1) source(1)],[tpoint(2) source(2)],[tpoint(3) source(3)],'Color',[0 0.003*count 0],'LineStyle','-','LineWidth',2)
%         tpoint = source;
        %------------------------------------------------------------------%
        T = T+t; % sum the total time
        if T > Time
            break
        end
        Dist = Dist + dist; % sum the total distance
        count = count + 1;
    end
end
toc
%% impulse response
out = sortrows(out,2);
Sample = Time*Fs;
TimePoints = 0:Sample-1;
IR = zeros(Sample,1);

for n = 1:1:N
    if out(n,2) ~= 0 && round(out(n,2)*Fs)<=Sample
%           IR = IR + amplitude{1,1}{1,n} * sinc(TimePoints-time{1,1}{1,n}*Fs).';
        t = TimePoints - round(out(n,2)*Fs);
        t(t==0) = Fs*2;
        t(t<Fs*2) = 0;
        t(t==max(t))= 1;
        IR = IR + out(n,1) * t.';
    end
end
% plot(TimePoints/Fs,IR,'LineWidth',2);
% title('Indoor Ray Tracing Impulse Response')
% xlabel('Time (s)')
% ylabel('Power')
% grid

%% plot the IR graph and compare it with real IR
if overlap_only == 1 || overlap_only == 0
    [x,fs] = audioread(diplay_audio);
    [ma,R]=max(x); % Find the direct soudn power, which is the lasrgest one
    R2 = find(IR~=0); % Find the first cell that has value
    di = R-R2(1); % find the position difference between the model and real IR
    figure;
    if overlap_only == 1
        % match the direct sound position
        if di > 0
            t = 1:Sample;
            x = x(di:Sample+di);
            plot((t)/Fs,x(t))
            hold on
            plot(TimePoints/Fs,IR,'LineWidth',2)
            title('Indoor Ray Tracing Impulse Response')
            xlabel('Time (s)')
            ylabel('Power')
            grid;
            hold off
        else
            t = 1:Sample;
            x = [zeros(-di,1);x(1:Sample+di)];
            plot((t)/Fs,x(t))
            hold on
            plot(TimePoints/Fs,IR,'LineWidth',2)
            title('Indoor Ray Tracing Impulse Response')
            xlabel('Time (s)')
            ylabel('Power')
            grid;
            hold off
        end
    % plot three graphs including overlapped echograph, real IR and model
    % IR
    else if overlap_only == 0
            % match the direct sound position
            if di > 0
                subplot(3,1,1)
                t = 1:Sample;
                x = x(di:Sample+di);
                plot((t)/Fs,x(t))
                hold on
                plot(TimePoints/Fs,IR,'LineWidth',2)
                title('Indoor Ray Tracing Impulse Response')
                xlabel('Time (s)')
                ylabel('Power')
                grid;
                hold off
            else
                subplot(3,1,1)
                t = 1:Sample;
                x = [zeros(-di,1);x(1:Sample+di)];
                plot((t)/Fs,x(t))
                hold on
                plot(TimePoints/Fs,IR,'LineWidth',2)
                title('Indoor Ray Tracing Impulse Response')
                xlabel('Time (s)')
                ylabel('Power')
                grid;
                hold off
            end
            subplot(3,1,2)
            plot((t)/Fs,x(t))
            subplot(3,1,3)
            plot(TimePoints/Fs,IR,'LineWidth',2)
            title('Indoor Ray Tracing Impulse Response')
            xlabel('Time (s)')
            ylabel('Power')
            grid;
        end
    end
end

% plot the room frame
% xx = zeros(wnum*size(wall,2),2);
% yy = zeros(wnum*size(wall,2),2);
% zz = zeros(wnum*size(wall,2),2);
% for u = 1:1:wnum % wall number counter
%     for v = 1:1:size(wall,2) % vertice counter for every wall
%         if (v == size(wall,2) && wall(u,v)~=0)||(wall(u,v) ~= 0 && wall(u,v+1)==0)
%             xx((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),1),vertex(wall(u,1),1)];
%             yy((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),2),vertex(wall(u,1),2)];
%             zz((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),3),vertex(wall(u,1),3)];
%             continue
%         else if v < size(wall,2) && wall(u,v+1)~=0
%                 xx((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),1),vertex(wall(u,v+1),1)];
%                 yy((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),2),vertex(wall(u,v+1),2)];
%                 zz((u-1)*size(wall,2)+v,:) = [vertex(wall(u,v),3),vertex(wall(u,v+1),3)];
%             end
%         end
%     end
% end
% xx = xx';
% yy = yy';
% zz = zz';
% plot3(xx, yy, zz,'k-','LineWidth',2)
% hold on
% % plot the source and the receiver position
% plot3(src(1),src(2),src(3),'.','MarkerSize',50) % Source position
% plot3(rcv(1),rcv(2),rcv(3),'.','MarkerSize',50) % Receiver position

