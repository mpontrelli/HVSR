%Path indicates the directory on a personal machine to where the codes have
%been placed.
%The code will process microtremor data on a specified station using the
%guidelines from Yilar et al. 2017

% Set the data time window and station
time1={'2011-11-23 17:00:00'}; %Start time of the event
time2={'2011-11-23 17:10:25'}; %End time of the event
Station='FFD'; %Stations and networks will be added
Network='NE';
Component={'BHZ','BHN','BHE'};

%Extracting data from the IRIS Network and performing analysis on each of
%the data
for ii = 1:length(time1)
  %Get data from IRIS
  [~,trace1UD,xUD] = getDMCData(char(time1{ii}),char(time2{ii}),Station,Network,Component{1});
  [~,~,xNS] = getDMCData(char(time1{ii}),char(time2{ii}),Station,Network,Component{2});
  [~,~,xEW] = getDMCData(char(time1{ii}),char(time2{ii}),Station,Network,Component{3});
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
     k(1) = k(2)+(40+25)*40;
  end 
  
  %Compute the complex time series
  xHmatrix = xNSmatrix + 1i.*xEWmatrix; 
  
  %Compute the fft for each 40 sec data window
  for iii = 1:10
      xUDmatrix(iii,:) = abs(fft(xUDmatrix(iii,:)))/1600;
      xNSmatrix(iii,:) = abs(fft(xNSmatrix(iii,:)))/1600;
      xEWmatrix(iii,:) = abs(fft(xEWmatrix(iii,:)))/1600;
      xHmatrix(iii,:) = abs(fft(xHmatrix(iii,:)))/1600;
  end
  
  %Stack each of the components
  [xUDavg, sigmaUD, confinthighUD, confintlowUD] =  HVSRavg(xUDmatrix);
  [xNSavg, sigmaNS, confinthighNS, confintlowNS] =  HVSRavg(xNSmatrix);
  [xEWavg, sigmaEW, confinthighEW, confintlowEW] =  HVSRavg(xEWmatrix);
  [xHavg, sigmaH, confinthighH, confintlowH] =  HVSRavg(xHmatrix);
  
  %Computing the x-axis
  N = length(xNSavg);
  fax_binsN = (0 : N-1); %samples in NS component
  fax_HzN1 = fax_binsN*fs/N; %frequency axis NS (Hz)
  if fs == 40
      ff=2;
  end
  N_2 = ceil(N/ff); %half magnitude spectrum
  fax_HzN = fax_HzN1(1 : N_2);
  
  %Smooth the data
  xUDavg = smooth(xUDavg,20); xNSavg = smooth(xNSavg,20); 
  xEWavg = smooth(xEWavg,20); xHavg = smooth(xHavg,20);
  
  %Compute the H/V ratio with the NS & EW component
  HV_NS = HV(xNSavg,xUDavg);
  HV_EW = HV(xEWavg,xUDavg);
  HV_H = HV(xHavg,xUDavg);
  
  %Plot the data
  figure 
  semilogy(fax_HzN,HV_EW(1:N_2))
  title('EW component')
  xlabel('Frequency (Hz)')
  ylabel('Amplification')
  
  figure 
  semilogy(fax_HzN,HV_NS(1:N_2))
  title('NS component')
  xlabel('Frequency (Hz)')
  ylabel('Amplification')
  
  figure 
  semilogy(fax_HzN,HV_H(1:N_2))
  title('Complex Time series')
  xlabel('Frequency (Hz)')
  ylabel('Amplification') 
end