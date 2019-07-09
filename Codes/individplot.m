function individplot(HV_final_matrix, newfaxhz, statname)
figure
plot(newfaxhz, HV_final_matrix, 'Color',  'k' , 'Linewidth', .5);
title(strcat(statname), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'FontSize',20,'YScale', 'log')
xlim([0 5])
ylim([0.1 100])
grid on    
end
