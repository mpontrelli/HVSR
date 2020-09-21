%% Hermes
% Gui for pulling seismic data from IRIS and looking at it 

%% Author: Marshall Pontrelli
% Date: 9/18/2020

% Here's a website to help you out
%https://www.mathworks.com/help/matlab/creating_guis/about-the-simple-programmatic-gui-example.html

function simple_gui2
% SIMPLE_GUI2 Select a data set from the pop-up menu, then
% click one of the plot-type push buttons. Clicking the button
% plots the selected data in the axes.

    %  Create and then hide the UI as it is being constructed.
    f = figure('Visible','off','Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
   
    %  Construct the components.
    htext = uicontrol('Style','text','String','Station','Fontsize', 12,'FontName', 'Times New Roman',...
        'Position',[10,945,100,30]);
    hpopup = uicontrol('Style','popupmenu',...
        'String',{'NE WES','CO HODGE', 'CO JSC'},'Fontsize', 12,'FontName', 'Times New Roman',...
        'Position',[10,925,100,25]);
    htext2 = uicontrol('Style','text','String','Frequency','Fontsize', 12,'FontName', 'Times New Roman',...
        'Position',[10,875,100,30]);
    hpopup2 = uicontrol('Style','popupmenu',...
        'String',{'Broadband','Short-Period', 'Long-Period'},'Fontsize', 12,'FontName', 'Times New Roman',...
        'Position',[10,850,100,25],...
        'Callback',@filterbutton_Callback);
    htext3 = uicontrol('Style','text','String','Channel','Fontsize', 12,'FontName', 'Times New Roman',...
        'Position',[10,800,100,30]);
    hpopup3 = uicontrol('Style','popupmenu',...
        'String',{'NS','EW', 'V'},'Fontsize', 12,'FontName', 'Times New Roman',...
        'Position',[10,775,100,25],...
        'Callback',@componentbutton_Callback);
    ha = axes('Units','Pixels','Position',[300,100,850,850]); 
   
    align([htext, hpopup],'Center','None');
    align([htext2, hpopup2],'Center','None');
    align([htext3, hpopup3],'Center','None');
    % Make the UI visible.
    % Initialize the UI.
    % Change units to normalized so components resize automatically.
    f.Units = 'normalized';
    ha.Units = 'normalized';
    hsurf.Units = 'normalized';
    hmesh.Units = 'normalized';
    hcontour.Units = 'normalized';
    htext.Units = 'normalized';
    hpopup.Units = 'normalized';
    htext2.Units = 'normalized';
    hpopup2.Units = 'normalized';
    htext3.Units = 'normalized';
    hpopup3.Units = 'normalized';

    % Pull Weston data
    datapath = strcat('C:\Users\',getenv('username'),'\Desktop\Data_pull\');
    load(strcat(datapath,'WES.mat'),'V_broad','NS_broad','EW_broad',...
        'NS_broad','EW_broad', 'V_broad',...
        'NS_short','EW_short','V_short','NS_long','EW_long','V_long','samplerate','sensitivity', 'sensunits');
    length_hour = samplerate*3600;
    num_win = ceil(length(V_broad)/length_hour);
    time_vec = linspace(0,60,length_hour);
    color_vec = [0 0 0; 1 0 0; 0 0 1; 0 0.4 0;0 0 0; 1 0 0; 0 0 1; 0 0.4 0;...
        0 0 0; 1 0 0; 0 0 1; 0 0.4 0;0 0 0; 1 0 0; 0 0 1; 0 0.4 0;0 0 0; 1 0 0; 0 0 1; 0 0.4 0;...
        0 0 0; 1 0 0; 0 0 1; 0 0.4 0;0 0 0; 1 0 0; 0 0 1; 0 0.4 0]; 
   
    %% now start plotting the waveforms
    for q = 1:num_win - 1
        color = color_vec(q,:);
        st = (q-1)*length_hour+1;
        se = q*length_hour;
        horz_pos = 1500 + (q-1)*2000;
        a = NS_broad(st:se)+horz_pos;
        hold on
        plot(time_vec, a, 'color', color);
        hold on
    end
    % plot last window
    color = color_vec(num_win,:);
    st = (num_win-1)*length_hour+1;
    horz_pos = 1500 + (num_win-1)*2000;
    a = NS_broad(st:end) + horz_pos;
    end_t = length(NS_broad) - ((num_win-1)*length_hour+1);
    samps = end_t/3600;
    time_end = linspace(0,samps,end_t +1);
    hold on
    plot(time_end, a, 'color', color);
    hold on
    xlim([0 60])
    ylim([0 24*2060])
    xticks([0 5 10 15 20 25 30 35 40 45 50 55 60])
    xticklabels({'0','5','10','15','20','25','30','35','40','45','50','55','60'})
    yticks([1500 3500 5500 7500 9500 11500 13500 15500 17500 19500 21500 23500,...
        25500 27500 29500 31500 33500 35500 37500 39500 41500 43500 45500 47500])
    yticklabels({'00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00',...
        '08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00',...
        '17:00','18:00','19:00','20:00','21:00','22:00','23:00'})
    ylabel('(UTC)','Fontsize', 18,'FontName', 'Times New Roman','Color','k')
    xlabel('Time (Minutes)', 'Fontsize', 18,'FontName', 'Times New Roman','Color','k')
    set(gca, 'YDir','reverse','Fontsize', 12,'FontName', 'Times New Roman')
    grid on
    hold off
    % Assign the a name to appear in the window title.
    f.Name = 'Simple GUI';

    % Move the window to the center of the screen.
    %movegui(f,'center')

    % Make the window visible.
    f.Visible = 'on';

    %  Pop-up menu callback. Read the pop-up menu Value property to
    %  determine which item is currently displayed and make it the
    %  current data. This callback automatically has access to 
    %  current_data because this function is nested at a lower level.
    function componentbutton_Callback(source,eventdata) 
        % Determine the selected data set.
        str = get(source, 'String');
        val = get(source,'Value');
%         r = get(handles.hpopup2, 'Value');
%         disp(r)
        switch str{val}
        case 'EW' % User selects EW.
            cla
        % Display plot.
            for q = 1:num_win - 1
                color = color_vec(q,:);
                st = (q-1)*length_hour+1;
                se = q*length_hour;
                 horz_pos = 1500 + (q-1)*2000;
                a = EW_broad(st:se)+horz_pos;
                hold on
                plot(time_vec, a, 'color', color);
                hold on
              end
            % plot last window
            color = color_vec(num_win,:);
            st = (num_win-1)*length_hour+1;
            horz_pos = 1500 + (num_win-1)*2000;
            a = EW_broad(st:end) + horz_pos;
            end_t = length(EW_broad) - ((num_win-1)*length_hour+1);
            samps = end_t/3600;
            time_end = linspace(0,samps,end_t +1);
            hold on
            plot(time_end, a, 'color', color);
            hold on
            xlim([0 60])
            ylim([0 24*2060])
            xticks([0 5 10 15 20 25 30 35 40 45 50 55 60])
            xticklabels({'0','5','10','15','20','25','30','35','40','45','50','55','60'})
            yticks([1500 3500 5500 7500 9500 11500 13500 15500 17500 19500 21500 23500,...
                25500 27500 29500 31500 33500 35500 37500 39500 41500 43500 45500 47500])
            yticklabels({'00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00',...
                '08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00',...
                '17:00','18:00','19:00','20:00','21:00','22:00','23:00'})
            ylabel('(UTC)','Fontsize', 18,'FontName', 'Times New Roman','Color','k')
            xlabel('Time (Minutes)', 'Fontsize', 18,'FontName', 'Times New Roman','Color','k')
            set(gca, 'YDir','reverse','Fontsize', 12,'FontName', 'Times New Roman')
            grid on
            hold off
        case 'NS' % User selects NS.
            cla
            % Display plot.
            for q = 1:num_win - 1
                color = color_vec(q,:);
                st = (q-1)*length_hour+1;
                se = q*length_hour;
                horz_pos = 1500 + (q-1)*2000;
                a = NS_broad(st:se)+horz_pos;
                hold on
                plot(time_vec, a, 'color', color);
                hold on
            end
            % plot last window
            color = color_vec(num_win,:);
            st = (num_win-1)*length_hour+1;
            horz_pos = 1500 + (num_win-1)*2000;
            a = NS_broad(st:end) + horz_pos;
            end_t = length(NS_broad) - ((num_win-1)*length_hour+1);
            samps = end_t/3600;
            time_end = linspace(0,samps,end_t +1);
            hold on
            plot(time_end, a, 'color', color);
            hold on
            xlim([0 60])
            ylim([0 24*2060])
            xticks([0 5 10 15 20 25 30 35 40 45 50 55 60])
            xticklabels({'0','5','10','15','20','25','30','35','40','45','50','55','60'})
            yticks([1500 3500 5500 7500 9500 11500 13500 15500 17500 19500 21500 23500,...
                25500 27500 29500 31500 33500 35500 37500 39500 41500 43500 45500 47500])
            yticklabels({'00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00',...
                '08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00',...
                '17:00','18:00','19:00','20:00','21:00','22:00','23:00'})
            ylabel('(UTC)','Fontsize', 18,'FontName', 'Times New Roman','Color','k')
            xlabel('Time (Minutes)', 'Fontsize', 18,'FontName', 'Times New Roman','Color','k')
            set(gca, 'YDir','reverse','Fontsize', 12,'FontName', 'Times New Roman')
            grid on
            hold off
       case 'V' % User selects V.
           cla
           % Display plot.
            for q = 1:num_win - 1
                color = color_vec(q,:);
                st = (q-1)*length_hour+1;
                se = q*length_hour;
                horz_pos = 1500 + (q-1)*2000;
                a = V_broad(st:se)+horz_pos;
                hold on
                plot(time_vec, a, 'color', color);
                hold on
            end
            % plot last window
            color = color_vec(num_win,:);
            st = (num_win-1)*length_hour+1;
            horz_pos = 1500 + (num_win-1)*2000;
            a = V_broad(st:end) + horz_pos;
            end_t = length(V_broad) - ((num_win-1)*length_hour+1);
            samps = end_t/3600;
            time_end = linspace(0,samps,end_t +1);
            hold on
            plot(time_end, a, 'color', color);
            hold on 
            xlim([0 60])
            ylim([0 24*2060])
            xticks([0 5 10 15 20 25 30 35 40 45 50 55 60])
            xticklabels({'0','5','10','15','20','25','30','35','40','45','50','55','60'})
            yticks([1500 3500 5500 7500 9500 11500 13500 15500 17500 19500 21500 23500,...
                25500 27500 29500 31500 33500 35500 37500 39500 41500 43500 45500 47500])
            yticklabels({'00:00','01:00','02:00','03:00','04:00','05:00','06:00','07:00',...
                '08:00','09:00','10:00','11:00','12:00','13:00','14:00','15:00','16:00',...
                '17:00','18:00','19:00','20:00','21:00','22:00','23:00'})
            ylabel('(UTC)','Fontsize', 18,'FontName', 'Times New Roman','Color','k')
            xlabel('Time (Minutes)', 'Fontsize', 18,'FontName', 'Times New Roman','Color','k')
            set(gca, 'YDir','reverse','Fontsize', 12,'FontName', 'Times New Roman')
            grid on
            hold off
        end % End Switch
    end % End component callback switch
end % End GUI function