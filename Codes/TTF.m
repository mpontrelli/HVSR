%% TTF
% compute the Vs30 of a velocity profile

    % INPUTS
    % d - vector of depths for soil profile
    
    % v - vector of velocities (same length as d) for soil profile
    
    % rho - vector of densities for soil profile
    
    % Q - vector of quality factors for soil profile (in units of 1/Q). Is equal to
    % 1/(2*zeta) where zeta is the damping ratio.
    
    % v_base - velocity for the basement (one number)
    
    % rho_base - density for the basement (one number)
    
    % num_freq - number of frequencies to compute
    
    % high_freq - highest frequency to be computed
    
    % depth_compute - the depth at which the response is computed (0 is the
    % free surface)
    
    % layer_half - layer of the halfspace, a big number means source is in
    % halfspcae
    
    % incidence - incidence angle (0 for vertically propagating S-waves)
    
    % OUTPUTS
    
    % vs_30 - Vs_30 of the profile
    
    % class - NEHRP site classification of the profile
    
function [freqs,amps] =  TTF(d,v,rho,Q, v_base,rho_base, num_freq,...
    high_freq,depth_compute,layer_half,incidence)
    % change into NRATTLE directory in HVSR code package
    cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes\NRATTLE'));
    
    % Now rewrite the ctl file with new profile
    [fid3,~] = fopen('nrattle.ctl','w');
    % Now write the file
    
    fprintf(fid3,'!Control file for nrattle (data lines are free format)\r\n');
    fprintf(fid3,'!Revision of program involving a change in the control file on this date:\r\n');
    fprintf(fid3,'   09/09/11\r\n');
    fprintf(fid3,'!Stem name for output files (the file containing model information and complex spectra will be named \r\n');
    fprintf(fid3,'! <stem>.nrattle_complex_spectra.out \r\n');
    fprintf(fid3,'! and the one containing the absolute spectral amplitudes used in plotting will be named \r\n');
    fprintf(fid3,'! <stem>.nrattle_amps4plot.out\r\n');
    fprintf(fid3,'!  nrattle_bj97_vs30_1100_ai_45\r\n');
    fprintf(fid3,'  test_nrattle_02mar11\r\n');
    fprintf(fid3,'!Number of frequencies (including f = 0.0) and highest freq. (nfreq, fhigh):   \r\n');
    a = strcat({'   '},num2str(num_freq),{'    '},num2str(high_freq));
    a = a{1};
    fprintf(fid3,strcat(a,'\r\n'));
    fprintf(fid3,'!Depth at which response is computed (0.0 for free surface):   \r\n');
    a = strcat({'    '},num2str(depth_compute));
    a = a{1};
    fprintf(fid3,strcat(a,'\r\n'));
    fprintf(fid3,'!Model information:\r\n');
    fprintf(fid3,'!For each layer- thickness, Vs, density, 1/Q (if > 1.0, the program assumes\r\n');
    fprintf(fid3,'!that the value represents Q; note that damping in percent = 100*0.5*1/Q):\r\n');
    % now loop through the model and get all the layer values
    for i = 1:length(d)
        qq = num2str(i);
        d_str = num2str(d(i));
        v_str = num2str(v(i));
        rho_str = num2str(rho(i));
        Q_str = num2str(Q(i));
        a = strcat({'    '},qq,{'   '},d_str,{'   '},v_str,{'   '},rho_str,{'   '},Q_str);
        a = a{1};
        fprintf(fid3,strcat(a,'\r\n'));
    end
    fprintf(fid3,'!  Halfspace Vs and density (1/Q is automatically set to 0.0):\r\n');
    rho_base_str = num2str(rho_base);
    v_base_str = num2str(v_base);
    a = strcat({'    '},v_base_str,{'   '},rho_base_str);
    a = a{1};
    fprintf(fid3,strcat(a,'\r\n'));
    fprintf(fid3,'!Layer corresponding to halfspace (can be less than actual halfspace)   \r\n');
    fprintf(fid3,'! (a large number means the source is in the halfspace) and \r\n');
    fprintf(fid3,'! angle of incidence (< 0.0 means obtain theta from site_amp file):    \r\n');
    layer_half_str = num2str(layer_half);
    incidence_str = num2str(incidence);
    a = strcat({'   '},layer_half_str,{'  '},incidence_str);
    a = a{1};
    fprintf(fid3,strcat(a,'\r\n'));
    %fprintf(fid3,f);
    fclose(fid3);
    % Run NRATTLE from matlab
    [~,~] = system('NRATTLE.exe < NRATTLE.ctl');
    
    % Now wait 3 seconds to let the files print
    pause(3)
    % now read the TTF
    filename='test_nrattle_02mar11.nrattle_amps4plot.out';
    M = dlmread(filename,'',19,0);
    freqs=M(:,1);
    amps=M(:,3);
    % change back to normal directory
    cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes'));
    figure
    plot(freqs,amps)
    title('Theoretical Transfer Function')
    xlabel('Frequency (Hz)','FontSize', 18)
    ylabel('Amplification','FontSize', 18)
    set(gca,'FontSize',20,'YScale', 'log')
    xlim([0 10])
    ylim([0.1 100])
    grid on
    set(gca,'YScale', 'log','FontName', 'Times New Roman', 'FontSize', 14)

    [maxamp,I]=max(amps);
    disp(strcat('maximum amplitude = ', num2str(maxamp)))
    FSF=freqs(I);
    disp(strcat('fundamental frequency = ', num2str(FSF)))
end