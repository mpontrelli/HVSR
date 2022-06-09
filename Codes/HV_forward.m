%% HV_forward

% Create a modeled HVSR curve from an S and P wave velocity profile and
% densities. This uses the engine from the HV-Inv software package
% developed by Antonio Garcia-Jerez and others

    % INPUTS
    
    % n - number of layers
    
    % t - thickenss vector (meters)

    % Vp - P-wave velocity vector
    
    % Vs - S-wave velocity vector

    % rho - density vector
    
    % Vp_base - P wave velocity of basement
    
    % Vs_base - S wave velocity of basement
    
    % rho_base - density of basement
    
    % fmin1 - low frequency value you want computed
    
    % fmax1 - high frequency value you want computed
    
    % NR - number of Rayleigh wave modes (10 is default in HV-Inv)
    
    % NL - number of Love wave modes (10 is default in HV-Inv)
    
    % dk - number of integration points (500 is default is HV-Inv)
    
    % log - 0 if you want linear spaced frequencies, 1 if you want log
    % spaced frequencies
    
    % OUTPUTS
    
    % HV - HVSR curve
    
    % f - Frequency vector
    
%% Author: Marshall Pontrelli
% Date: 10/29/2021

function [HV,f] = HV_forward(n,t, Vp, Vs, rho, Vp_base, Vs_base, rho_base, fmin1, fmax1, NR, NL, dk, log)

    %% input apsv
    apsv = 1e-3;
    
    %% Make general paths which allow HV-Inv to run
    %% Script principal for execute GUI HV-INV
    addpath(genpath('./etc'));
    addpath(genpath('./bin/Functions/Func_HV'));
    addpath(genpath('./bin/Functions/Func_Inv/F_GUI'));
    addpath(genpath('./bin/Functions/Func_Inv/F_Plot'));
    addpath(genpath('./bin/Functions/Func_Inv/F_Kernel/F_Fow'));
    addpath(genpath('./bin/Functions/Func_Inv/F_Kernel/F_Inv'));
    addpath(genpath('./bin/Functions/Func_Inv/F_Kernel/F_Mod'));
    %% create m matrix
    m = zeros(n+2, 4);
    m(1,1) = n;
    m(2:end - 1,1) = t;
    m(2:end - 1, 2) = Vp;
    m(2:end - 1, 3) = Vs;
    m(2:end - 1, 4) = rho;
    m(n+2,1) = 0;
    m(n+2,2) = Vp_base;
    m(n+2,3) = Vs_base;
    m(n+2,4) = rho_base;
    
    %% Now move to the folder with the HV forward model function from HV-Inv
    cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes\HV-INV-master\HV-Inv_source_MATLAB\bin\Functions\Func_INV\F_Kernel\F_Fow'));
    
    %% And compute the HVSR model
    [HV,f] = HVPD(m,fmin1,fmax1,n,NR,NL,dk,log,apsv);
    %% Now move back to the HVSR codes folder
    cd(strcat('C:\Users\',getenv('username'),'\Desktop\HVSR\Codes'));
end