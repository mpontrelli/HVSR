close all
clear all
%% EW1
filename = 'C:\Users\mpontr01\Box\Projects\Kik-Net\Data\TKCH08\TKCH080010070819.EW2';

[fid3,~] = fopen(filename,'r');

% Now read the file
%
% Skip the first 10 lines 
for jjj = 1:10
    line = fgetl(fid3);
end

% read fs
fs = line(19:21);

for jjj = 1:4
    line = fgetl(fid3);
end
b = length(line);
scale_factor_num = str2num(line(19:22));
scale_factor_denom = str2num(line(29:b));
scale_factor = scale_factor_num/scale_factor_denom;

%Set space deliminatorIt
deliminator='';
%Data starts in row 109, column 1
R=17; %Row data start text file
C=0; %column data start text file
data=dlmread(filename,deliminator,R,C); %retrieve data
data = data';
xNS=data(:);
xNS = scale_factor*(xNS - mean(xNS));
figure
plot(xNS)
% xV=data(:,2); %Vertical_Component
% xV = scale_factor*(xV - mean(xV));
% figure
% plot(xV)
% xEW=data(:,3); %East_West_Component
% xEW = scale_factor*(xEW - mean(xEW));
% figure
% plot(xEW)
% 
% Y = data(:,4); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% Y = data(:,5); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% Y = data(:,6); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)


% %% EW2
% filename = 'C:\Users\mpontr01\Box\Projects\Kik-Net\Data\TKCH08\TKCH080010070819.EW2';
% 
% [fid3,~] = fopen(filename,'r');
% 
% % Now read the file
% %
% % Skip the first 10 lines 
% for jjj = 1:10
%     line = fgetl(fid3);
% end
% 
% % read fs
% fs = line(19:21);
% 
% for jjj = 1:4
%     line = fgetl(fid3);
% end
% b = length(line);
% scale_factor_num = str2num(line(19:22));
% scale_factor_denom = str2num(line(29:b));
% scale_factor = scale_factor_num/scale_factor_denom;
% 
% %Set space deliminatorIt
% deliminator='';
% %Data starts in row 109, column 1
% R=17; %Row data start text file
% C=0; %column data start text file
% data=dlmread(filename,deliminator,R,C); %retrieve data
% xNS=data(:,1); %North_South_Component
% xNS = scale_factor*(xNS - mean(xNS));
% figure
% plot(xNS)
% xV=data(:,2); %Vertical_Component
% xV = scale_factor*(xV - mean(xV));
% figure
% plot(xV)
% xEW=data(:,3); %East_West_Component
% xEW = scale_factor*(xEW - mean(xEW));
% figure
% plot(xEW)
% 
% Y = data(:,4); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% Y = data(:,5); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% Y = data(:,6); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% %% NS1
% filename = 'C:\Users\mpontr01\Box\Projects\Kik-Net\Data\TKCH08\TKCH080010070819.NS1';
% 
% [fid3,~] = fopen(filename,'r');
% 
% % Now read the file
% %
% % Skip the first 10 lines 
% for jjj = 1:10
%     line = fgetl(fid3);
% end
% 
% % read fs
% fs = line(19:21);
% 
% for jjj = 1:4
%     line = fgetl(fid3);
% end
% b = length(line);
% scale_factor_num = str2num(line(19:22));
% scale_factor_denom = str2num(line(29:b));
% scale_factor = scale_factor_num/scale_factor_denom;
% 
% %Set space deliminatorIt
% deliminator='';
% %Data starts in row 109, column 1
% R=17; %Row data start text file
% C=0; %column data start text file
% data=dlmread(filename,deliminator,R,C); %retrieve data
% xNS=data(:,1); %North_South_Component
% xNS = scale_factor*(xNS - mean(xNS));
% figure
% plot(xNS)
% xV=data(:,2); %Vertical_Component
% xV = scale_factor*(xV - mean(xV));
% figure
% plot(xV)
% xEW=data(:,3); %East_West_Component
% xEW = scale_factor*(xEW - mean(xEW));
% figure
% plot(xEW)
% 
% Y = data(:,4); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% Y = data(:,5); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% Y = data(:,6); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% %% NS2
% filename = 'C:\Users\mpontr01\Box\Projects\Kik-Net\Data\TKCH08\TKCH080010070819.NS2';
% 
% [fid3,~] = fopen(filename,'r');
% 
% % Now read the file
% %
% % Skip the first 10 lines 
% for jjj = 1:10
%     line = fgetl(fid3);
% end
% 
% % read fs
% fs = line(19:21);
% 
% for jjj = 1:4
%     line = fgetl(fid3);
% end
% b = length(line);
% scale_factor_num = str2num(line(19:22));
% scale_factor_denom = str2num(line(29:b));
% scale_factor = scale_factor_num/scale_factor_denom;
% 
% %Set space deliminatorIt
% deliminator='';
% %Data starts in row 109, column 1
% R=17; %Row data start text file
% C=0; %column data start text file
% data=dlmread(filename,deliminator,R,C); %retrieve data
% xNS=data(:,1); %North_South_Component
% xNS = scale_factor*(xNS - mean(xNS));
% figure
% plot(xNS)
% xV=data(:,2); %Vertical_Component
% xV = scale_factor*(xV - mean(xV));
% figure
% plot(xV)
% xEW=data(:,3); %East_West_Component
% xEW = scale_factor*(xEW - mean(xEW));
% figure
% plot(xEW)
% 
% Y = data(:,4); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% Y = data(:,5); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% Y = data(:,6); %East_West_Component
% Y = scale_factor*(Y - mean(Y));
% figure
% plot(Y)
% 
% 
