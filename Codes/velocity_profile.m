%% velocity_profile
% plot a shear wave velocity profile from input vectors of depth and 
% velocity

    % INPUTS
    % h - vector of layer thicknesses
    
    % v - vector of shear wave velocities correlating to the thicknesses

%% Author: Marshall Pontrelli
% Date: 4/20/2020


function velocity_profile(h, v)
    h_new = 0;
    for i = 1: length(h) % Sum up the thicknesses
        h_new(i+1) = h_new(i) + h(i);
    end
    h_new = repelem(h_new,2);
    h_new(1) = []; 
    h_new(end) = [];
    v_new = repelem(v,2); 
    plot(v_new, h_new, 'linewidth', 2, 'color','k')
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 18, 'Ydir','reverse');
    xlabel('Velocity')
    ylabel('Depth')
    xlim([0 v_new(end) + v_new(1)])
    grid on
    box on
end