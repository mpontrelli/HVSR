%% Plot all SSRs Mexico_city
close all
clear all


sitelist = {'AE02','AL01','AO24', 'AP68','AR14','AU11','AU46','BA49','BL45','BO39',...
    'CA20', 'CA59', 'CB43','CE18', 'CE23','CE32','CH84','CI05','CJ03',...
    'CJ04','CO47', 'CO56','CP28','CT64', 'CU80', 'DM12','DR16', 'DX37','EO30','ES57','EX08','EX09','EX12','GA62',...
    'GC38','GR27','HA41','HJ72', 'IB22','JA43','JC54','LI33', 'LI58', 'LV17','ME52',...
    'MI15', 'MY19', 'NZ20', 'NZ31','PA34', 'PD42','PE10', 'RI76', 'RM48',...
    'SI53', 'SP51', 'TH35', 'TL08', 'TL55', 'UC44', 'VG09', 'VM29',...
    'XP06'};

reflist = {'CS78', 'FJ74', 'IM40','TE07', 'TP13', 'UI21'};

fold_path = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\SSR_Shape_statistics\';

% Load frequency vector
load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\AE02\AE0219900511234349')
fax_HzN = data.processing.filtereddata.freq_vec;
upbound = 10;
lowbound = 0.1;
[~, lowbound] = min(abs(fax_HzN - lowbound));
[~, upbound] = min(abs(fax_HzN - upbound));

for ii = 1:length(sitelist)
    statname = sitelist{ii};
    fig = figure;
    xlim([0.1 10])
    ylim([0.1 100])
    xticks([.1 1 10])
    xticklabels({'0.1', '1', '10'})
    yticks([1 10 100])
    yticklabels({ '1','10', '100'})
    xlabel('Frequency (Hz)')
    ylabel('Amplification')
    title(statname)
    set(gca,'YScale', 'log','XScale','log', 'FontName', 'Times New Roman', 'FontSize', 18)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    grid on
    box on
    hold on
    
    % CS78
    load(strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\SSR_Shape_statistics\',...
    statname,'\',statname,'-CS78'))
    ahatfSSR = data.complex.ahatf';
    CS78 = plot(fax_HzN, ahatfSSR, 'Color', 'm');
    hold on
    
    % FJ74
    load(strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\SSR_Shape_statistics\',...
    statname,'\',statname,'-FJ74'))
    ahatfSSR = data.complex.ahatf';
    FJ74 = plot(fax_HzN, ahatfSSR, 'Color', 'c');
    hold on
    
    % IM40
    load(strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\SSR_Shape_statistics\',...
    statname,'\',statname,'-IM40'))
    ahatfSSR = data.complex.ahatf';
    IM40 = plot(fax_HzN, ahatfSSR, 'Color', 'r');
    hold on
    
    % TE07
    load(strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\SSR_Shape_statistics\',...
    statname,'\',statname,'-TE07'))
    ahatfSSR = data.complex.ahatf';
    TE07 = plot(fax_HzN, ahatfSSR, 'Color', 'g');
    hold on
    
    % TP13
    load(strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\SSR_Shape_statistics\',...
    statname,'\',statname,'-TP13'))
    ahatfSSR = data.complex.ahatf';
    TP13 = plot(fax_HzN, ahatfSSR, 'Color', 'b');
    hold on
    
    % UI21
    load(strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\SSR_Shape_statistics\',...
    statname,'\',statname,'-UI21'))
    ahatfSSR = data.complex.ahatf';
    UI21 = plot(fax_HzN, ahatfSSR, 'Color', 'k');
    hold on
    
    % HVSR
    load(strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Shape_statistics\',statname));
    ahatfHV = shapedata.complex.ahatf';
    HV = plot(fax_HzN, ahatfHV, 'Color', 'k', 'LineStyle','--', 'Linewidth', 2);
    hold on

    legend([CS78, FJ74, IM40,TE07, TP13, UI21, HV],'CS78', 'FJ74', 'IM40','TE07', 'TP13', 'UI21', 'HVSR')
    hold off
    saveas(fig, strcat('C:\Users\mpontr01\Box\Geohazards Research Group\2020_5_6\All_SSRs\',statname,'.jpg'))
    close all
end