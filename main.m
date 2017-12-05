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

% 125 250 500 1000 2000 4000
% 1   2   3    4    5    6
Band_num = 4;

% Initial_MMR
% Initial_Pollack
Initial_Tanna
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
% imprange = [4]; % set the second order impulse that we want to see 0==all -1==none
refl_order = 2; % reflection order 1,2,-1

N = 270*270;            % Set up the number of the rays emitted must be power of 2
% r = log10(V)*norm(src-rcv)*sqrt(4/N);    % Set up the radius of listener
r = 0.4;
delta = 0.005;        % Set up the absoption coefficient of the walls
alpha = 0.001;        % Set up the air absorption coefficient
c = 331.4 + 0.6*23.6; % Set up the sound speed: c + 0.6T, where T is room temperature
Po = 1;               % Set up the power of the sound
P_thresh = 0.001*Po;  % Set up the threshhold of the ray energy
Time = 0.1;               % Time 

% generate the function of every wall
plane = TPlane(wall,wnum,vertex);

% Here I initialize a matrix storing the direction of each ray
ray_direc = ray_direction(N);
% Set up output storing matrix
out = zeros(N,2+refl_order);
        
%% Reflection computation
for n = 1:1:N
    source = src;
    newdirec = ray_direc(n,:);
    P = Po;
    Dist = 0;
    T = 0;
    count = 0;
    while(count <= refl_order)
        % see if the ray is reaching the listener
        temp = rcv-source;
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
            out(n,1) = P* exp(-alpha*dtemp2)*spk_coef*2*dtemp1/(4/3*pi*r^3);
            out(n,2) = Ts;
            break
        end
        [source,newdirec,dist,P,t,out] = refl(source,newdirec,plane,wall,vertex,P,alpha,c,out,n,count,beta);
        % [source,newdirec,dist,P,t] = refl(source,newdirec,plane,wall,vertex,P,beta,alpha,c,f_concrete,f_wood,f_lime);
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
% out(1,:) = 0;
% out(1,1) = E_direc;
% out(1,2) = T_direc;
% out = sortrows(out,2);
dt = norm(src-rcv)/c;
Sample = ceil((Time+dt)*fs);
TimePoints = 0:Sample-1;
IR = zeros(Sample,1);

for n = 1:1:N
    if out(n,2) ~= 0 && round(out(n,2)*fs)<=Sample
%           IR = IR + amplitude{1,1}{1,n} * sinc(TimePoints-time{1,1}{1,n}*Fs).';
        t = TimePoints - round(out(n,2)*fs);
        t(t==0) = fs*2;
        t(t<fs*2) = 0;
%         t(t==max(t))= 1;
        t_start = find(t==max(t));
        temp = out(n,1);
        IR(t_start:t_start+size(temp,2)-1) = IR(t_start:t_start+size(temp,2)-1)+ temp;
    end
end
%% plot the IR graph and compare it with real IR
if overlap_only == 1 || overlap_only == 0
    
    y = filter(B,A,y,[],1);
    IR = filter(B,A,IR,[],1);
    
    [ma,R]=max(y); % Find the direct soudn power, which is the lasrgest one
    IR=IR*ma/max(IR); % normalize the model energy and the measured energy
    R2 = find(IR~=0); % Find the first cell that has value
    di = R-R2(1); % find the position difference between the model and real IR
    
    figure;
    if overlap_only == 1
        % match the direct sound position
        if di > 0
            t = 1:Sample;
            y = y(di:Sample+di);
            plot((t)/fs,y(t))
            hold on
            plot(TimePoints/fs,IR(1:Sample),'LineWidth',1)
            title('Indoor Ray Tracing Impulse Response')
            xlabel('Time (s)')
            ylabel('Power')
            grid;
            hold off
        else
            t = 1:Sample;
            y = [zeros(1,-di) y(1:Sample+di)];
            plot((t)/fs,y(t))
            hold on
            plot(TimePoints/fs,IR,'LineWidth',1)
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
                y = y(di:Sample+di);
                plot((t)/fs,y(t))
                hold on
                plot(TimePoints/fs,IR,'LineWidth',1)
                title('Indoor Ray Tracing Impulse Response')
                xlabel('Time (s)')
                ylabel('Power')
                grid;
                hold off
            else
                subplot(3,1,1)
                t = 1:Sample;
                y = [zeros(-di,1);y(1:Sample+di)];
                plot((t)/fs,y(t))
                hold on
                plot(TimePoints/fs,IR,'LineWidth',1)
                title('Indoor Ray Tracing Impulse Response')
                xlabel('Time (s)')
                ylabel('Power')
                grid;
                hold off
            end
            subplot(3,1,2)
            plot((t)/fs,y(t))
            subplot(3,1,3)
            plot(TimePoints/fs,IR,'LineWidth',1)
            title('Indoor Ray Tracing Impulse Response')
            xlabel('Time (s)')
            ylabel('Power')
            grid;
        end
    end
end



