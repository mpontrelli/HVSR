close all
clear all

depth = [-4.96 -5.83 -8.60 -6.2832 -5.99];
dist = [0 50 325 400 500];
vert1x = [0,0];
vert1y = [0, max(depth)];
vert2x = [max(dist), max(dist)];
vert2y = [0, -5.99];
surfx = [0, max(dist)];
surfy = [0, 0];
basinplot = figure;
hold on
fill([dist flip(dist)],[depth zeros(size(depth))], [0.772549 0.35294 0.0666667],'LineStyle','none')
plot(dist,depth, 'Color', 'k', 'LineWidth', 2)
plot(vert1x, vert1y, 'Color', 'k', 'LineWidth', 2)
plot(vert2x,vert2y, 'Color', 'k', 'LineWidth', 2)
plot(surfx, surfy, 'Color', 'k', 'LineWidth', 2)


plot(dist(1), 0, 'o','MarkerSize', 15, 'MarkerEdgeColor', [0.5176 0 0.6588], 'MarkerFaceColor', [0.5176 0 0.6588])
% 
% %Lake Stations
plot(dist, zeros(length(dist)), 'o','MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');

%makes figure full screen
set(gcf, 'Units', 'Normalized', 'OuterPosition', [1, 0.5, 1, 0.5]);
ylabel('Depth (m)')
xlabel('distance (m)')

set(gca,'FontSize',20)
set(gca,'fontname','times')
annotation('textbox',[.13 .82 .1 .2],'String','1a','EdgeColor','none', 'FontSize', 20, 'fontname', 'times')
annotation('textbox',[.22 .82 .1 .2],'String','2a','EdgeColor','none', 'FontSize', 20, 'fontname', 'times')
annotation('textbox',[.6 .82 .1 .2],'String','6a','EdgeColor','none', 'FontSize', 20, 'fontname', 'times')
annotation('textbox',[.71 .82 .1 .2],'String','7a','EdgeColor','none', 'FontSize', 20, 'fontname', 'times')
annotation('textbox',[.87 .82 .1 .2],'String','8a','EdgeColor','none', 'FontSize', 20, 'fontname', 'times')
saveas(basinplot, 'C:\Users\mpontr01\Desktop\boston_site_response\Presentations\10_31_2019\basinplot.jpg');

