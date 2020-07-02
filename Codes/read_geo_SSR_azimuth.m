%% Plot HVSR Mexico_city
close all
clear all

sitelist = { 'AE02','AL01','AO24', 'AP68','AR14','AU11','AU46','BA49','BL45','BO39',...
     'CA20', 'CA59', 'CB43','CC55', 'CE23','CE32','CH84','CI05','CJ03',...
     'CJ04','CO47', 'CO56', 'CU80', 'DM12','DR16', 'DX37','EO30','ES57','EX08','EX09','EX12','GA62',...
     'GC38','GR27','HJ72', 'IB22','JA43','JC54','LI33', 'LI58', 'LV17','ME52',...
    'MI15', 'MY19', 'NZ20', 'NZ31', 'PD42','PE10', 'RI76', 'RM48',...
    'SI53', 'SP51', 'TH35', 'TL08', 'TL55', 'UC44', 'VG09', 'VM29',...
    'XP06'};
% ref = 'FJ74';
sheet = 8;
outpath = strcat('C:\Users\mpontr01\Box\Pontrelli_et_al_2020\SSR_W\');
load('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2\AE02\AE0219900511234349')
filename = 'C:\Users\mpontr01\Box\Pontrelli_et_al_2020\Final_tables.xlsx';
fax_HzN = data.processing.filtereddata.freq_vec;
upbound = 10;
lowbound = 0.1;
[~, lowbound] = min(abs(fax_HzN - lowbound));
[~, upbound] = min(abs(fax_HzN - upbound));
for ii = 1:length(sitelist)
    statname = sitelist{ii};
    %if strcmp(statname,ref) ==0 %&& strcmp(statname,'CC55') ==0
        rownum = num2str(ii + 1);
        load(strcat('C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\geomean_SSR_W\',...
        statname))

        ahatf = data.ahatf';
        sigma = data.sigma';
        confinthigh = data.confinthigh';
        confintlow = data.confintlow';
        number_events = data.num;

        %% Now quantify the shape of the peaks
        [peakfreqs, peakamps, hpb, f1s, f2s, Areamat, proms, amps, peakind2, freqs, sigs, I1s, I2s] = peakiden(ahatf, fax_HzN, sigma, lowbound, upbound);
        [~, I] = max(amps);

        %% Now set some conditions to classify the peak
        [~, freq_class] = sort(freqs, 'descend');
        [~, amp_class] = sort(amps, 'descend');
        [~, prom_class] = sort(proms, 'descend');
        [~, hpb_class] = sort(hpb, 'descend');
        [~, area_class] = sort(Areamat, 'descend');
        [~, sig_class] = sort(sigs, 'descend');
        classmatrix = vertcat(freq_class, amp_class, prom_class, hpb_class, area_class, sig_class);

        %% Now save peak data
        datamat = vertcat(freqs(I),amps(I),proms(I),hpb(I),Areamat(I),sigs(I), number_events);
        datamat = datamat';
        %% now plot, make the figure and set the base
        fin_fig = figure;
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

        %% start with the confidence interval
        fr = fax_HzN(lowbound:length(fax_HzN))';
        cohr = confinthigh(lowbound:length(fax_HzN))';
        colr = confintlow(lowbound:length(fax_HzN))';
        confidenceinterval=shadedplot(fr, cohr, colr,[.9,.9,.9],[1 1 1]);
        hold on
    
        %% Plot filled in peaks color coordinated based on their area
        freq_ar = peakfreqs{I};
        amp_ar = peakamps{I};

        fill(freq_ar, amp_ar, 'b', 'LineStyle','none','FaceAlpha',0.5)
        hold on



                freq_ar = peakfreqs{I};
                amp_ar = peakamps{I};
                ampd = amp_ar(I1s(I));
                ampdd = amp_ar(I2s(I));
                plot([f1s(I),f2s(I)],[ampd, ampdd], 'c', 'LineWidth', 2)
                hold on



        %% now plot fundamental resonance

                line([freqs(I),freqs(I)],[0.1, 100],'LineStyle', '--', 'color','k')
                hold on




        %% Now the hpb-prominence cross
        hold on

                plot([freqs(I),freqs(I)],[amps(I), amps(I) - proms(I)], 'c', 'Linewidth', 2)
                hold on


        %% Now plot ahatf
        plot(fax_HzN,ahatf, 'LineWidth', 2, 'Color', [0 0.5 0])
        %% Now print some important info
        if length(amps) > 0
            if length(amps) ==1
                a = "Max peak freq =" + " "  + num2str(freqs,3);
                b = strcat("Amp =" + " "  + num2str(amps,3));
                c = strcat("HPB =" + " "  +num2str(hpb,3));
                d = strcat("Prom =" + " "  +num2str(proms,3));
                e = strcat("Area =" + " "  +num2str(Areamat,3));
                f = strcat("\sigma =" + " "  +num2str(sigs,3));
                str = {a, b, c, d, e, f};
                ylim=get(gca,'ylim');
                xlim=get(gca,'xlim');
                q = find(ahatf == amps);
                if fax_HzN(q) < 1
                    text(xlim(2)-5,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
                else
                    text(xlim(1)+0.3,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
                end
                else
                    [~, I] = max(amps);
                    a = "Max peak freq =" + " "  + num2str(freqs(I),3);
                    b = strcat("Amp =" + " "  + num2str(amps(I),3));
                    c = strcat("HPB =" + " "  +num2str(hpb(I),3));
                d = strcat("Prom =" + " "  +num2str(proms(I),3));
                e = strcat("Area =" + " "  +num2str(Areamat(I),3));
                f = strcat("\sigma =" + " "  +num2str(sigs(I),3));
                str = {a, b, c, d, e, f};
                ylim=get(gca,'ylim');
                xlim=get(gca,'xlim');
                q = find(ahatf == amps(I));
                if fax_HzN(q) < 1
                    text(xlim(2)-5,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
                else
                    text(xlim(1)+0.3,ylim(2)-60,str, 'FontName', 'Times New Roman', 'FontSize', 18,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
                end

            end
        else
            str = 'No peaks';
            text(0.3,15,str, 'FontName', 'Times New Roman', 'FontSize', 24,'Color', 'black', 'HorizontalAlignment', 'right', 'EdgeColor','k','BackgroundColor', 'w')
    
        end
% 
         saveas(fin_fig, strcat(outpath,statname,'.jpg'))
%         datamat_fin = datamat(:,1)';
% 
%     
        writematrix(datamat,filename,'Sheet',sheet,'Range',strcat('D',rownum,':J',rownum))
        clear xlim
        clear ylim
        close all
    end


% str = {strcat('Peak 1 prom = ',{' '},num2str(proms(1))), strcat('Peak 2 prom = ',{' '},num2str(proms(2)))};
% stra = str{1}{1};
% strb = str{2}{1};
% strc = strcat('Peak 1 area = ',num2str(Areamat(1)));
% strd = strcat('Peak 1 hpb = ',num2str(hpbss(1)));
% stre = strcat('Peak 2 area = ',num2str(Areamat(2)));
% strf = strcat('Peak 2 hpb = ',num2str(hpbss(2)));
% strreal = {stra, strd, strc};
% text(.15,20,strreal, 'FontName', 'Times New Roman', 'FontSize', 12,'Color', 'red')
% strreal2 ={strb, strf,stre};
% text(.15,5,strreal2, 'FontName', 'Times New Roman', 'FontSize', 12, 'Color', 'blue')
% 
% % Now the arrow for freq and amp
% freq = strcat('fn = ',num2str(matrix(2,1)));
% amp = strcat('Amp = ',num2str(matrix(2,2)));
% strreal2 = {freq, amp};
% text(freqs(2)-5 ,amps(2) +1,strreal2,'FontName', 'Times New Roman', 'FontSize', 12, 'Color', 'blue' )
