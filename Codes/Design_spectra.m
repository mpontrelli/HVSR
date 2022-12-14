%% Design_spectra

    % Compute a design spectra. Use this document to see how it's done https://www.graitec.com/Help/Advance-Design/desktop/Response_spectrum_curves_(ASCE).pdf
    % 

    % INPUTS
    % Ss - Parameter from USGS Maximum Credible Earthquake map. the mapped MCER spectral response acceleration parameter at short periods as
           %determined in accordance with Section 11.4.1

    % S1 - Parameter from USGS Maximum Credible Earthquake map. the mapped MCER spectral response acceleration parameter at a period of 1 s as
           %determined in accordance with Section 11.4.1

    % Fa - Short period site coefficient

    % Fv - Long period site coefficient

    % Tl - Long period transition period

    % MCE site: https://earthquake.usgs.gov/hazards/designmaps/pdfs/?code=IBC&edition=2012
    
    
  %% Author: Marshall Pontrelli
  % Written Summer 2020 for my qualifying exam, updated and turned into a
  % function Summer 2022 while writing my dissertation.

%% Start  
function Design_spectra(Ss,S1,Fa,Fv,Tl)
    Sms = Fa*Ss;
    Sm1 = Fv*S1;
    Sds = (2/3)*Sms;
    Sd1 = (2/3)*Sm1;
    T0 = 0.2*(Sd1/Sds);
    Ts = Sd1/Sds;

    % loop over period
    for ii = 1:1000
        Tn = (ii-1)*0.01;
        if Tn < T0
            Sd(ii) = Sds*(0.4+0.6*Tn / T0);
        end
        if Tn >= T0 && Tn <= Ts
            Sd(ii) = Sds;
        end
        if Tn > Ts && Tn <= Tl
            Sd(ii) = Sd1/Tn;
        end
        if Tn > Tl 
            Sd(ii) = Sd1*Tl/Tn^2;
        end
        period(ii) = Tn;
    end
   
    %% now plot
    figure
    plot(period,Sd, 'linewidth',2);
    xlim([0 5])
    grid on
    box on
    xlabel('Period (secs)')
    ylabel('Spectral acceleration (g)')
    title('5% damped response spectra')
    set(gca,'FontName', 'Times New Roman', 'FontSize', 18)
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
end