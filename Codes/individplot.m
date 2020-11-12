function individplot(HV_final_matrix, newfaxhz, statname, sav,outpath)
individualfiltered = figure;
plot(newfaxhz, HV_final_matrix, 'Color',  'k' , 'Linewidth', .5);
title(strcat(statname), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
xlim([0 20])
ylim([0.1 100])
grid on   

if strcmp(sav, 'yes') == 1
    saveas(individualfiltered, strcat(outpath, '\', 'individualfiltered.jpg'));
end
end
