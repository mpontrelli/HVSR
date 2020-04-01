%% plot vecotrs of sigma
close all
clear all
% Author Marshall Pontrelli
% Date: 11/20/2019

codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Shape_statistics';

cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
sigma_mat = [];
freq_mat = [];
for i = 1: length(stationlist)
    station = stationlist(i);
    filename = strcat(station.folder, '\', station.name);
    load(filename)
    sigma_mat(i,:) = shapedata.complex.sig_vec;
    shapes2 = shapedata.complex.shapes;
    freq = cell2mat(shapes2(2,2));
    freq_mat(i,:) = freq;
end
[B, I] = sort(freq_mat);
sigma_freq = shapedata.sigma_freq;
cd(codepath);

%% now smooth
for iii = 1:10
    sigma_mat2(iii,:) = smooth(sigma_mat(iii,:),400);
end

sigma_mat3 = sigma_mat2(I,:);
%% now plot
figure

for i = 1:length(stationlist)
    color = 1/(2*i);
    color = [color color color];
    plot(sigma_freq, sigma_mat3(i,:), 'Color',  color , 'Linewidth', .5);
    hold on
end
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('\sigma','FontSize', 18)
set(gca,'XScale', 'log', 'FontName', 'Times New Roman', 'FontSize', 14)
xlim([0.1 10])
xticks([.1 1 10])
xticklabels({'0.1', '1', '10'})
ylim([0.1 1])
title('Mexico City  RACM station \sigma values')
grid on   
