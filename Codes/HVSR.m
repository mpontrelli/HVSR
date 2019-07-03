function HVSR(path, datapath)
%close all
%close all
warning('off','all') %The warnings are from the triangular filter which is 
%still a piece of the code, though it can be removed. 

%% AU11
 stationlist = {'AU11'};%,'AL01', 'AO24', 'AP68', 'AR14', 'AU11', 'AU46', 'BA49',...
%     'BL45', 'BO39', 'CA20', 'CA59', 'CB43', 'CC55',...
%     'CE18', 'CE23', 'CE32', 'CH84', 'CI05', 'CJ03', 'CJ04', 'CO47', 'CO56',...
%     'CP28', 'CS78', 'CT64', 'CU80', 'DM12', 'DR16', 'DX37',...
%     'EO30', 'ES57', 'EX08', 'EX09', 'EX12', 'FJ74', 'GA62', 'GC38', 'GR27',...
%     'HA41', 'HJ72', 'IB22', 'IM40', 'JA43', 'JC54', 'LI33', 'LI58', 'LV17',...
%     'ME52', 'MI15', 'MY19', 'NZ20', 'NZ31', 'PA34', 'PD42', 'PE10', 'RI76',...
%     'RM48', 'SI53', 'SP51', 'TE07', 'TH35', 'TL08',...
%     'TL55', 'UC44', 'VG09', 'VM29', 'XP06'};
for eee = 1:length(stationlist)
    station = stationlist{eee};
    disp(station)
    d = strcat(datapath,'\',station);
    %go into data directory and build structure of all files in it
    cd(d)
    %cd 'C:\Users\Marshall\Box Sync\tFall_2018\Research\Mexico_City\Data\AE02';
    files = dir;
    files = files(3:length(files));
    
    %create empty matrix that gets filled with all H/V values, counter is used
    %to index this matrix
    HV_final_matrix = [];
    peakfreq = [];
    peakamp = [];
    counter = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %function

    %change directory back to codes to access functions needed 
    %cd 'C:\Users\Marshall\Box Sync\tFall_2018\Research\Mexico_City\Codes';
    d = strcat(path, '\HVSR\Codes');
    cd(d)
    for file = files'
        filename = strcat(datapath,'\', station,'\',file.name);
        %filename = strcat('C:\Users\Marshall\Box Sync\tFall_2018\Research\Mexico_City\Data\AE02\',file.name);
        [xNS,xV,xEW, fs] = readfile1(filename);
        [PGANS,PGAV,PGAEW] = PGA(xNS,xV,xEW);
        if PGANS < 0.1 && PGAV < 0.1 && PGAEW < 0.1
            counter = counter + 1;
            %fs = station.fs; %sampling frequency in hz
            [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
            [N_2, fax_HzN, XH_magfilt,XV_magfilt] =  Magresp(xNS, xV, xEW, fs); %Compute mag responses and run through triangular filter
    
            %perform H/V
            [H_V1] = HV(XH_magfilt,XV_magfilt);

            %make Hz vector and linear interpolate all H/V ETFs to this vector
            newfaxhz = 0:0.001:20;
            newH_V1 = interp1(fax_HzN, H_V1, newfaxhz);
            HV_final_matrix(counter, :) = newH_V1; 
            clear H_V1
            else
            [xNS, xV, xEW] = Butter(xNS, xV, xEW, fs); %filter the data
            [N_2, fax_HzN, XH_magfilt,XV_magfilt] =  Magresp(xNS, xV, xEW, fs); %Compute mag responses and run through triangular filter
    
            %perform H/V
            [H_V1] = HV(XH_magfilt,XV_magfilt);
   
            %make Hz vector and linear interpolate all H/V ETFs to this vector
            newfaxhz = 0:0.001:20;
            newH_V1 = interp1(fax_HzN, H_V1, newfaxhz);
%
%     figure
%     plot(newfaxhz, newH_V1, 'Color', [0 .30196 .6588], 'LineWidth', 1.5)
%     title(file.name)
%     xlabel('frequency (Hz)')
%     ylabel('HVSR')
%     set(gca, 'YScale', 'log')
%     grid on
%     txt = {strcat('PGA NS = ',num2str(PGANS)),strcat('PGA EW = ', num2str(PGAEW))};
%     text(1.25,30,txt,'FontSize',16)
%
            clear newfaxhz
            clear newH_V1
            clear H_V1 
        end
    end
newfaxhz = 0:0.001:20;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%statistics per Thompson et al 2012 page 34
%compute maximum likelihood estimator of median

[ahatf, sigma, confinthigh, confintlow] = HVSRavg(HV_final_matrix);

HVSRplot(ahatf, newfaxhz, sigma, confinthigh, confintlow, station);

%Call nrattle and plot TTF
%cd 'C:\Users\mpontr01\Box Sync\Box Sync\tFall_2018\Research\Mexico_City\Data\CE32_data';

% filename = 'TTF_workbook1.xlsx';
% % sheet = 1;
% % xlRange = 'E2:GJM2';
% TTF = xlsread(filename,1,'E7: GJM7');
% TTF=TTF(10:length(TTF)-1);
% 
% TTFplot = plot(newfaxhz,TTF, 'Linewidth', 1.5, 'Color', 'k');
%legend([a],{'SH1D TTF'})
%saveas(fig,strcat(output,'CE32test.pdf'))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find median sigmali between 1st and 4th peaks



%%
%Now we run through all the significant peaks and compute shape statistics
% 
% Taxstat = [];
% for f = 1:length(peakamp)
%     taxstat(f,1) = f; 
%     taxstat(f,2) = peakfreq(f); 
%     A = peakamp(f);
%     taxstat(f,3) = A;
%     amploc2 = amplocs2(f);
%     [I1, I2, f1, f2, hpb] =  HalfPowerBand2(A, amploc2, newfaxhz, ahatf); 
% %     disp(f1)
% %     disp(f2)
%     taxstat(f,4) = hpb; 
%     a = sigma(I1:I2);
%     sigmai = median(a);
%     taxstat(f,5) = sigmai;
% end
%points = plot(freqpeak,amps,'o','MarkerSize',7,'MarkerEdgeColor','g','MarkerFaceColor','g');

% b = ahatf((locs1(1)):(locs1(4)));
% c = TTF((locs1(1)):(locs1(4)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%copute sigmai and pearson's correlation coefficient between first and
%fourth peaks
% sigmai = median(a);
% rrr = num2str(sigmai);
% rrrr = rrr(1:5);
% [r,p] = corrcoef(b,c);
% sss = num2str(r(2));
% ssss = sss(1:5);
% if sigmai < 0.35 && r(2) > 0.60
%     cat = 'LG';
% end
% if sigmai < 0.35 && r(2) < 0.60
%     cat = 'LP';
% end
% if sigmai > 0.35 && r(2) < 0.60
%     cat = 'HP';
% end
% if sigmai > 0.35 && r(2) > 0.60
%     cat = 'LG';
% end
% txt = {strcat('sigma = ',rrrr),strcat('r = ',ssss), 'Depth = 74.98 m', 'Damping = 1.7%'};
% text(1.25,30,txt,'FontSize',16)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%test peaks figure
% figure
% set(gca, 'YScale', 'log')
% hold on
% plot(b)
% plot(c)

%saveas(ETF,strcat(output,'\',station,'.jpeg')) 
N = length(ahatf); %length North_South_Component
width = .1; %width for triangle moving average filter in hz
q = ceil((N/20)*width); %width for triangle moving average filter in samples
%e = smoothdata(ahatf, 'gaussian', q);
% ee = smooth(e,q);
% eef = smooth(ee,q);
% eeg = smooth(eef,q);
% eeh = smooth(eeg,q);
% eei = smooth(eeh,q);
% eej = smooth(eei,q);
% eek = smooth(eej,q);
% ETF = plot(newfaxhz, ahatf, 'Color', [0 0.30196 .6588] , 'Linewidth', 3);
%plot(newfaxhz,e, 'LineWidth', 2, 'Color', 'r')
%plot(newfaxhz,ee, 'LineWidth', 2, 'Color', 'g')
%plot(newfaxhz,eek, 'LineWidth', 2, 'Color', 'k')

%%
% Determine if peak is a peak
count = 0;
[amps,amplocs] = findpeaks(e);
for hh = 1:length(amps)
    freqpeak(hh) = newfaxhz(amplocs(hh));
end
% plot(freqpeak, amps, 'o', 'color', 'k')
% if amps(1) < 3
%     amps = amps(2:length(amps));
% end
% if length(locs1)>=4
[trough,troughloc] = findpeaks(-1*e);
for hh = 1:length(trough)
    freqtrough(hh) = newfaxhz(troughloc(hh));
end
if freqtrough(1) > freqpeak(1)
    Y = 1;
    X = 1;
    trough = horzcat(Y,trough);
    troughloc = horzcat(X,troughloc);
end

if freqpeak(length(freqpeak)) > freqtrough(length(freqtrough))
    Y = -1;
    X = 19991;
    trough = horzcat(trough, Y);
    troughloc = horzcat(troughloc,X);
end
for hh = 1:length(trough)
    freqtrough(hh) = newfaxhz(troughloc(hh));
end
% plot(freqtrough, -1 * trough, '*', 'color', 'r')
for k = 1:length(amps)
    height = amps(k)/sqrt(2);
    lefttrough = -1 * trough(k);
    righttrough = -1 * trough(k + 1);
%     disp(num2str(amps(k)))
%     disp(num2str(lefttrough))
%     disp(num2str(righttrough))
    if lefttrough < height && righttrough < height
        count  = count + 1;
        peakfreq(count) = newfaxhz(amplocs(k));
        peakamp(count) = amps(k); 
        amplocs2(count) = amplocs(k);
    end
end

points = plot(peakfreq,peakamp,'o','MarkerSize',15,'MarkerEdgeColor','g','MarkerFaceColor','g');
%% Compute desired statistics

Taxstat = [];
for f = 1:length(peakamp)
    taxstat(f,1) = f; 
    taxstat(f,2) = peakfreq(f); 
    A = peakamp(f);
    taxstat(f,3) = A;
    amploc2 = amplocs2(f);
    [I1, I2, f1, f2, hpb] =  HalfPowerBand2(A, amploc2, newfaxhz, ahatf); 
%     disp(f1)
%     disp(f2)
    taxstat(f,4) = hpb; git p
    a = sigma(I1:I2);
    sigmai = median(a);
    taxstat(f,5) = sigmai;
end
hold off

fclose('all')
end
end