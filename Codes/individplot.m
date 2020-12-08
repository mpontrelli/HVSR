function individplot(HV_final_matrix, newfaxhz, statname, sav,lowbound,upbound,outpath)
individualfiltered = figure;
plot(newfaxhz, HV_final_matrix, 'Color',  'k' , 'Linewidth', .5);
title(strcat(statname), 'FontSize', 20)
xlabel('Frequency (Hz)','FontSize', 18)
ylabel('Amplification','FontSize', 18)
set(gca,'YScale', 'log', 'XScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)
ylim([0.1 100])
grid on 
xlim([newfaxhz(lowbound) newfaxhz(upbound)])
%xlim([fax_HzN(1) 40])
ylim([0.1 100])
xticks([0.1 1 10])
xticklabels({'0.1','1','10', num2str(upbound)})
yticks([0.1 1 10 100])
yticklabels({'0.1','1','10', '100'})

if strcmp(sav, 'yes') == 1
    saveas(individualfiltered, strcat(outpath, '\', 'individualfiltered.jpg'));
end
end
