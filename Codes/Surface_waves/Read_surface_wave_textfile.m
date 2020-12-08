%% Read_surface_wave_textfile

% Reads and plots the text file output from the RAS-24 system to be used
% for surface wave analysis

%% Author: Marshall Pontrelli
% Date: 9/11/2020
close all
clear all

%% inputs
min_offset = 2;
max_offset = 24;
fs = 500;
datapath = 'C:\Users\mpontr01\Box\Projects\Surface_wave\9_14\';
%% import the data
% ACCESSING THE DATA
% go into the data folder and get a list of stations
path = pwd;
cd(datapath)
tracelist = dir;
tracelist = tracelist(3:length(tracelist));
% start the for loop that goes through all the station folders
B = zeros(501, 24);
for eee = 1:length(tracelist)
    trace = tracelist(eee); % read the folder info
    tracename = trace.name;
    filename = strcat(datapath,tracename);
    A = importdata(filename);
    B = A + B;
    q{eee} = importdata(filename);
end
B = B/length(tracelist);

%% figure out how big B is

B = B(:,1:23);
c = size(B);
a = c(2);
b = c(1);
time = 1000*(0:b-1)/(fs);

%% Now normalize B
for i = 1:a
    sta = B(:,i);
    max_sta = max(sta);
    col = sta/max_sta;
    B_norm(:,i) = col;
end


%% now plot
figure
for i = 1:a
    sta = B_norm(:,i);
    % Extract positive and negative part
    yp = (sta + abs(sta))/2;
    yn = (sta - abs(sta))/2;
    hold on
    plot(time,yn+i,'b')
    hold on
    fill(time,yp+i,i,'Facecolor','k')
    hold on
end
view([90, 90])
grid on
box on
xlabel('Time (ms)','FontName', 'Times New Roman', 'FontSize', 18,'rotation', 270, 'VerticalAlignment','middle');

yticks([1 a])
yticklabels({num2str(min_offset),num2str(max_offset)})
ax = gca;
ax.YAxisLocation = 'right';


set(gca, 'FontName', 'Times New Roman', 'FontSize', 18);    
xlim([0,500])
ylim([-0.5, i+1])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 0.5, 0.96]);
title('Offset (m)','FontName', 'Times New Roman', 'FontSize', 18);

%% Now fft
figure
for i = 1:a
    sta = B_norm(:,i);
    X = abs(fft(sta));
    X_phase = unwrap(angle(fft(sta)));
    %Computing the frequency -axis
    N = length(X);
    fax_binsN = (0 : N-1); %samples in NS component
    fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
    N_2 = ceil(N/2); %half magnitude spectrum
    fax_HzN = fax_HzN1(1 : N_2);
    X_mag = X(1 : N_2);
    X_phase = X_phase(1 : N_2);
    
    hold on
    plot(fax_HzN,X_mag)
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);


%% now make a dispersion curve
% X = fft(B,b,1)/b; % fft along columns (time)
% P = angle(X);
% A = abs(X);
% norm = X./A;
% V = P.*norm;
% V_new = horzcat(V(:,1),zeros(b,a-1));
% for i = 1:a-1
%     col = V(:,i+1);
%     V_new(:,i) = V_new(:,i) + col;
% end
% Vf = abs(V_new);
% figure
% heatmap(Vf)


%% Edit the parameters below that represent your dataset:

%load SampleData.mat % Load the synthetic dataset
U=B_norm; % my x-t domain signals defined for this script
pick='manual'; % This sets the function for picking the dispersion curve. 
             % Pick setting options are 'manual' , 'auto' or 'none'
waterlevel=0; % Dispersion Image waterlevel expressed as a percent for display/picking 
              % purposes. Amplitudes in the dispersion image that are less
              % than "waterlevel" percent of the peak amplitudes will be 
              % set to NANs. E.g., a waterlevel=0 will allow all Dispersion
              % Image data to be displayed, while a waterlevel=95 will
              % only allow the top 5 percent of amplitudes to be
              % displayed. High values may be desired, especially when
              % 'manual' picking is set. 
            
ct=50:1:1200; % A vector for the velocity of interest in m/s.
               % This is the velocity "scanning" range
               
freqlimits=[5 80]; % For plotting purposes only. Sets the frequency (Hz)
                    % plotting range for the dispersion images. Enter ['none']
                    % if you desire the natural limits based on time sampling

% Input Signal Properties %
fs=500; % sampling frequency Hz
min_x = 2; % min offset [m] of reciever spread
max_x = 24; % max offset [m] of reciever spread
d_x = 1; % distance between receivers [m]

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% End Inputs and Run Script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = min_x:d_x:max_x; % spatial sampling vector (receiver position)
t=1/fs:1/fs:length(U(:,1))/fs;  % time vector in seconds

%%%%% Find the dispersion image %%%%%%
%% 1. First, find the fft of each trace
% determine the proper frequency vector
L=length(t); % number os samples
T=max(t); % max time
if ~mod(L,2)
    f = (1/T)*(0:L/2-1); % n is even
else
    f = (1/T)*(0:(L-1)/2); % n is odd
end

Ufft=fft(U); % fft is performed on each column (the time vector)
Ufft=Ufft(1:length(f),:); % only keep the fft results for the positve frequencies
w=f*2*pi; % convert the freq vector to rads/sec

%% 2. Normalize the spectrum
Rnorm=Ufft./abs(Ufft); % Normalize the spectrum

%% 3. Perfrom the summation over N in the frequency domain
for ii=1:length(w) % send in frequencies one at a time
    As=zeros(length(x),length(ct)); %initialize As or overwite the previous As
    for n=1:length(x) % send in a positon (x) one at a time
        % Iterate As over position
        As(n,:)=exp(1j*w(ii)*(x(n)./ct))*Rnorm(ii,n);
        % Note: The sign convention of MATLAB FFT is negative, so the minus
        % sign is dropped from the "-j" seen in Ryden 2004
    end
     AsSum(ii,:)=abs(sum(As)); % Sum up As for a given frequency. 
end                              % The final result will be a matrix of 
                                 % amplitudes changing with phase 
                                 % velocity in the x direction and 
                                 % frequency in the y direction

AsSum=AsSum'; % transpose the matrix so velocity is on vertical and freq is on horizontal
normed=AsSum./max(AsSum);% normalize the dispersion image so max column values=1

%%
%%%% This next section will pick the dispersion curve base on the "pick"
%%%%%%% setting ('auto' 'manual' or 'none')

%%%%%%%%%%%% curve autopicking %%%%%%%%%%%%%%%%%%%
% This section will auto pick the dispersion curve based on the peak 
% value in each column. This will need updating when 
% higher modes are introduced
if strcmp('auto',pick) == 1 % if auto picking is set
    remaindat=normed; % create a new variable called remaindat from the normalized image
    remaindat(remaindat < 1) = 0;% set all data amplitudes in the new variable that are less than 1 to zero
    [row,col]=find(remaindat); % find the indices where there is non-zero data (the peak that we set equal to 1)
    autovel=ct(row); % autopicked velocity indices that contain the peak value
    autofreq=f(col); % autopicked frequency indices that contain the peak value
    % get rid of all the extra picks at zero frequency
    nz=numel(autofreq)-nnz(autofreq); % determine the number of repeated picks at zero
    autofreq=autofreq(nz:end); % remove the repeated zero picks
    autovel=autovel(nz:end); % remove the repeated zero picks
    DispersionVelocity=autovel; % need to get the rlowess working again, I no longer have access to the toolbox that allows it
    normed(normed < waterlevel/100)=nan; % throw out data less than waterlevel
%%%%%%%%%%%%%%%% end autopicking %%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%% manual curve picking %%%%%%%%%%%%%
elseif strcmp('manual',pick) == 1 % if manual picking is set
    normed(normed < waterlevel/100)=nan; % throw out data less than waterlevel
    figure('units','normalized','outerposition',[0 0 1 1])
    h=pcolor(f,ct,normed);
    set(h, 'EdgeColor', 'none');
    colormap([jet;1 1 1]);
    ylabel('Phase Velocity (m/s)')
    xlabel('Freq (Hz)')
        if strcmp('none',freqlimits) == 1

        else
            xlim(freqlimits)
        end
    title({'1) Use your mouse to pick the dispersion curve using as many points as possible.';...
    '2) Press ENTER when finished ';'3) MATLAB will interpolate between points';'DUPLICATE PICKS WILL RESULT IN AN ERROR.'})
    set(gca,'FontSize',15,'fontweight','bold')
    [pick_f,pick_ct]=ginput;
    DispersionVelocity=interp1(pick_f,pick_ct,f,'spline');
    hold on
    plot(pick_f,pick_ct,'ko','markersize',10,'linewidth',5)
    plot(f,DispersionVelocity,'k','linewidth',2)
    title('Press ENTER to continue')
    pause
    close
%%%%%%%%% end manual curve picking %%%%%%%%%%%%

% No curve picking
elseif strcmp('none',pick) == 1 % if no picking is desired i.e., 'none' is set
    DispersionVelocity=nan(length(f));
    normed(normed < waterlevel/100)=nan; % throw out data less than waterlevel
end

DispersionCurve=[f' DispersionVelocity']; % this only has meaning if the curve was picked using 'manual' or 'auto'

%%
%%% End calculations. Plot everything %%%%
% Plot the normalized x-t data
figure;
subplot(1,2,1)
datanorm=U./max(U);
imagesc(x,t,datanorm)
colormap(jet)
xlabel('Receiver Position (m)')
set(gca,'xaxisLocation','top')
ylabel('time (sec)')
%title('x-t domain')
%ylim([0 0.05])
set(gca,'FontSize',25,'fontweight','bold')

% Plot the dispersion results
subplot(1,2,2)
h=pcolor(f,ct,normed); % plot the normalized dispersion image
set(h, 'EdgeColor', 'none');
colormap(jet)
hold on
plot(DispersionCurve(:,1),DispersionCurve(:,2),'k','linewidth',2) % Plot the picked dispersion curve. This is nothing if 'none' picking was set
if strcmp('none',freqlimits) == 1 % check for plotting limits
    % if 'none' then ignore the xlim
else
    xlim(freqlimits) % set the xlim if specified
end
ylabel('Phase Velocity (m/s)')
xlabel('Freq (Hz)')
title('Dispersion Image')
set(gca,'FontSize',25,'fontweight','bold')

