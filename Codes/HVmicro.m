
%Path indicates the directory on a personal machine to where the codes have
%been placed.
%The code will process microtremor data on a specified station using the
%guidelines from Yilar et al. 2017

%Extract the Station names with the network assocaited 
%cd(datapath)
%Table = readtable('Useable_Stations_NE.xls','Sheet','All_Stations');
%Network = Table(:,1); Station = Table(:,2); Station = table2array(Station);
%Network = table2array(Network);
%cd(path)

% Set the data time window and station
%time1={'2011-11-23 17:00:00'}; %Start time of the event
%time2={'2011-11-23 17:10:25'}; %End time of the event
time1={'2012-06-15 08:00:00'}; %Start time of the event
time2={'2012-06-15 08:10:25'}; %End time of the event
Station='FFD'; %Stations and networks will be added
Network='NE';
Component={'BHZ','BHN','BHE'};

%Extracting data from the IRIS Network and performing analysis on each of
%the data
for vi = 1:1
    for ii = 1:length(time1)
      %Get data from IRIS
      [~,trace1UD,xUD] = getDMCData(char(time1{ii}),char(time2{ii}),Station,Network,Component{1});
      [~,~,xNS] = getDMCData(char(time1{ii}),char(time2{ii}),Station,Network,Component{2});
      [~,~,xEW] = getDMCData(char(time1{ii}),char(time2{ii}),Station,Network,Component{3});
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
      k = [1,40];
      for iii = 1:10
         xUDmatrix(iii,:) = xUD((k(1)):(k(2)*40+k(1))-1);
         xNSmatrix(iii,:) = xNS((k(1)):(k(2)*40+k(1))-1);
         xEWmatrix(iii,:) = xEW((k(1)):(k(2)*40+k(1))-1);
         k(1) = k(1)+2600;
      end

      %Compute the complex time series
      xHmatrix = xNSmatrix + 1i.*xEWmatrix; 

      %Compute the fft for each 40 sec data window
      windowlen = 10;
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
      if fs == 40
          ff=2;
      end
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
      for iii = 1:10
          [H_VEW(iii,:)] = HV(XEWmatrix(iii,:),XUDmatrix(iii,:));
          [H_VNS(iii,:)] = HV(XNSmatrix(iii,:),XUDmatrix(iii,:));
          [H_VH(iii,:)] = HV(XHmatrix(iii,:),XUDmatrix(iii,:));
      end

      %Half magnitude spectrum
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
      title1 = strcat('EW Component',{' '},Station); 
      title2 = strcat('NS Component',{' '},Station);
      title3 = strcat('Complex Time Series',{' '},Station);
      HVSRplot(ahatfEW, fax_HzN, confinthighEW, confintlowEW, title1)
      [maxEW,positionEW] = max(ahatfEW); HzEW = fax_HzN(positionEW);
      %cd(record)
      %d = strcat(title1,'.pdf');
      %saveas(gcf,d);
      %cd(path)
      HVSRplot(ahatfNS, fax_HzN, confinthighNS, confintlowNS, title2)
      [maxNS,positionNS] = max(ahatfNS); HzNS = fax_HzN(positionNS);
      %cd(record)
      %d = strcat(title2,'.pdf');
      %saveas(gcf,d);
      %cd(path)
      HVSRplot(ahatfH, fax_HzN, confinthighH, confintlowH, title3)
      [maxH,positionH] = max(ahatfH); HzH = fax_HzN(positionH);
      %cd(record)
      %d = strcat(title3,'.pdf');
      %saveas(gcf,d);
      %cd(path)
      
      % Determine if peak is a peak
      count = 0;
      width = .5; %width for triangle moving average filter in hz
      window = ceil((N/20)*width);
      e = smooth(ahatfH,window);
      [amps,amplocs] = findpeaks(e);
      for hh = 1:length(amps)
            freqpeak(hh) = fax_HzN(amplocs(hh));
      end
      [trough,troughloc] = findpeaks(-1*ahatfH);
      for hh = 1:length(trough)
          freqtrough(hh) = fax_HzN(troughloc(hh));
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
      for hh = 1:length(troughloc)-1
          freqtrough(hh) = fax_HzN(troughloc(hh));
      end
      for k = 1:length(amps)
          height = amps(k)/sqrt(2);
          lefttrough = -1 * trough(k);
          righttrough = -1 * trough(k + 1);
          if lefttrough < height && righttrough < height
             count  = count + 1;
             peakfreq(count) = fax_HzN(amplocs(k));
             peakamp(count) = amps(k); 
             amplocs2(count) = amplocs(k);
          end
      end
      figure
      plot(fax_HzN, ahatfH)
      hold on
      plot(fax_HzN(amplocs2),peakamp,'go')

    end
end 