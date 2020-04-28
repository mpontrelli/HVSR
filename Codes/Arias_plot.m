%% Arias_plot 

    % INPUTS
    
    % OUTPUTS
    
    
%% Author: Marshall Pontrelli
% Date: 4/17/2020
% edited - 2/11/2020 - formatting and axis labeling

function Arias_plot(time, Ianorm)
    figure
    [~,I5] = min(abs((Ianorm - 0.05)));
    D5 = time(I5);
    [~,I95] = min(abs((Ianorm - 0.95)));
    D95 = time(I95);
    Brackdur = D95 - D5;
    [~,I90] = min(abs((Ianorm - 0.90)));
    FNplot = plot(time, Ianorm);
    sigdurFN=line([D5 D5], [0 max(Ianorm)],'color','g');
    line([D95 D95], [0 max(Ianorm)],'color','g')
    title('Normalized Arias Intensity for combined ground motions')
    xlabel('Time (secs)')
    ylabel('Normalized Arias Intensity')
    grid on
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14);
end