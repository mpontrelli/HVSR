record = 'C:\Users\justi\Desktop\Research\Record\USarrayHV';
path = 'C:\Users\justi\Desktop\Research\HVSR\Codes';
datapath = 'C:\Users\justi\Box\Justin Summer 2019\Preliminary Work\Station List';
database = 'Useable_Stations_NE.xls';
%Path indicates the directory on a personal machine to where the codes have
%been placed.
%The code will process microtremor data on a specified station using the
%guidelines from Yilar et al. 2017

%Extract the Station names with the network assocaited 
cd(datapath)
Table = readtable(database);
Network = Table(:,1); Station = Table(:,2); Station = table2array(Station);
Network = table2array(Network);
cd(path)

% Set the data time window and station
%time1={'2011-11-23 17:00:00'}; %Start time of the event
%time2={'2011-11-23 17:10:25'}; %End time of the event
time1={'2011-04-15 08:00:00'}; %Start time of the event
time2={'2011-04-15 08:10:25'}; %End time of the event
%Station = {'BRYW', 'EMMW', 'FFD', 'HNH', 'PQI', 'QUA2', 'WES', 'WVL', 'TRY'};
Component={'BHZ','BHN','BHE'};

%Extracting data from the IRIS Network and performing analysis on each of
%the data
for iv = 1:length(Station)
    for ii = 1:length(time1)
      %Get data from IRIS
      [~,trace1UD,xUD] = getDMCData(char(time1{ii}),char(time2{ii}),Station{iv},Network{iv},Component{1});
      [~,~,xNS] = getDMCData(char(time1{ii}),char(time2{ii}),Station{iv},Network{iv},Component{2});
      [~,~,xEW] = getDMCData(char(time1{ii}),char(time2{ii}),Station{iv},Network{iv},Component{3}); 
      %Conditional Statements to skip stations that do not have a record of
      %certian time requested
      if xUD(1) == 0 && xNS(1) == 0 && xEW(1) == 0
        continue
      end 
      %Stations that just recorded the vertical component
      if length(xUD) ~= length(xNS) || length(xUD) ~= length(xEW)
        continue
      end 
      fs = trace1UD.sampleRate;

      %Window the data with a non-overlapping window of 40 secs and 25 secs
      %apart, so collecting 625 secs of data in total
      k = [1,fs];
      xUDmatrix = []; xNSmatrix = []; xEWmatrix = [];
      for iii = 1:10
         xUDmatrix(iii,:) = xUD((k(1)):(k(2)*fs+k(1))-1);
         xNSmatrix(iii,:) = xNS((k(1)):(k(2)*fs+k(1))-1);
         xEWmatrix(iii,:) = xEW((k(1)):(k(2)*fs+k(1))-1);
         k(1) = k(1)+25*fs;
      end

      %Compute the complex time series
      xHmatrix = xNSmatrix + 1i.*xEWmatrix; 

      %Compute the fft for each 40 sec data window
      windowlen = 10; XUDmatrix = []; XEWmatrix = []; XNSmatrix = []; 
      XHmatrix = [];
      for iii = 1:windowlen
          XUDmatrix(iii,:) = abs(fft(xUDmatrix(iii,:)))/1600;
          XNSmatrix(iii,:) = abs(fft(xNSmatrix(iii,:)))/1600;
          XEWmatrix(iii,:) = abs(fft(xEWmatrix(iii,:)))/1600;
          XHmatrix(iii,:) = abs(fft(xHmatrix(iii,:)))/1600;
      end

      %Computing the x-axis
      [m,n] = size(xUDmatrix);
      N = n;
      fax_binsN = (0 : N-1); %samples in NS component
      fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
      ff = fs/20;
      N_2 = ceil(N/ff); %half magnitude spectrum
      fax_HzN = fax_HzN1(1 : N_2);

      %Smooth the magnitude Response
      N = n; %length HVSR
      width = .1; %width for triangle moving average filter in hz
      window = ceil((N/20)*width); %width for triangle moving average filter in samples where 20 is the number of Hz on your x-axis
      for iii = 1:windowlen
          XUDmatrix(iii,:) = smooth(XUDmatrix(iii,:),window);
          XNSmatrix(iii,:) = smooth(XNSmatrix(iii,:),window);
          XEWmatrix(iii,:) = smooth(XEWmatrix(iii,:),window);
          XHmatrix(iii,:) = smooth(XHmatrix(iii,:),window);
      end

      %Compute the HVSR
      H_VEW = []; H_VNS = []; H_VH = []; 
      for iii = 1:10
          [H_VEW(iii,:)] = HV(XEWmatrix(iii,:),XUDmatrix(iii,:));
          [H_VNS(iii,:)] = HV(XNSmatrix(iii,:),XUDmatrix(iii,:));
          [H_VH(iii,:)] = HV(XHmatrix(iii,:),XUDmatrix(iii,:));
      end

      %Half magnitude spectrum
      XH_VEW = []; XH_VNS = []; XH_VH = [];
      for iii = 1:windowlen
        XH_VEW(iii,:) = H_VEW(iii,1:N_2);
        XH_VNS(iii,:) = H_VNS(iii,1:N_2);
        XH_VH(iii,:) = H_VH(iii,1:N_2);
      end 

      %Average the HVSR
      [ahatfEW, sigmaEW, confinthighEW, confintlowEW] =  HVSRavg(XH_VEW);
      [ahatfNS, sigmaNS, confinthighNS, confintlowNS] =  HVSRavg(XH_VNS);
      [ahatfH, sigmaH, confinthighH, confintlowH] =  HVSRavg(XH_VH);
      
      %Plot the data & save the plots
      lowbound = 1; %Nyquist = fs/4;
      titleH = [Station{iv},' Complex'];
      HVSRplot(ahatfH, fax_HzN, confinthighH, confintlowH, lowbound, Station{iv})
      title(titleH)
      cd(record)
      saveas(gcf,titleH)
      cd(path)
      
      titleEW = [Station{iv},' East-West'];
      HVSRplot(ahatfEW, fax_HzN, confinthighEW, confintlowEW, lowbound, Station{iv})
      title(titleEW)
      cd(record)
      saveas(gcf,titleEW)
      cd(path)
      
      titleNS= [Station{iv},' North-South'];
      HVSRplot(ahatfEW, fax_HzN, confinthighEW, confintlowEW, lowbound, Station{iv})
      title(titleNS)
      cd(record)
      saveas(gcf,titleNS)
      cd(path)
      
      % Determine if peak is a peak
      [matrixH,peakind,~,~] = peakiden(ahatfH, fax_HzN, lowbound);
      %[taxstat] = specratstat(peakind, matrixH, ahatfH, fax_HzN, sigmaH, Station{iv}, lowbound);     
    end
end 