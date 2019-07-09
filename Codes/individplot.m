function individplot(HV_final_matrix, newfaxhz, station)
plot(newfaxhz(10:length(newfaxhz)-1), HV_final_matrix(10:length(HV_final_matrix)-1), 'Color',  'k' , 'Linewidth', .5);
title(strcat(station), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([0 5])
ylim([0.1 100])
grid on    
end
