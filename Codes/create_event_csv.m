%% create_event_csv
  % Reads the data structure from create_mat_file and turns it into a .csv
  % file
%% Author: Marshall Pontrelli
% Date: 4/16/2020


close all
clear all

%% ACCESS THE DATA
% go into the data folder and get a list of stations
codepath = 'C:\Users\mpontr01\Desktop\HVSR\Codes';
datapath = 'C:\Users\mpontr01\Box\Data\Ground motion\Mexico CIty\Processed_data2';

cd(datapath)
stationlist = dir;
stationlist = stationlist(3:length(stationlist));
counter = 1;

% Create matrix "A" and title the columns
A = {'event_ID', 'event_mag', 'event_lat', 'event_lon', ...
    'event_depth', 'event_azimuth', 'event_epi_distance',...
    'station_name', 'station_lat', 'station_lon', 'station_alt',...
    'station_soil', 'instrument_model', 'instrument_frequency',...
    'NS_PGA', 'A_NS_arias_intensity','A_NS_arias_D595', ...
    'A_NS_arias_D575', 'A_NS_arias_rate', 'A_NS_spectra_max',...
    'A_NS_spectra_max_period','A_NS_spectra_seconds02','A_NS_spectra_seconds05', ...
    'A_NS_spectra_seconds1', 'A_NS_spectra_seconds2', 'A_NS_spectra_seconds5',...
    'A_NS_HV_amp','A_NS_HV_freq',...
    'EW_PGA', 'A_EW_arias_intensity','A_EW_arias_D595', ...
    'A_EW_arias_D575', 'A_EW_arias_rate', 'A_EW_spectra_max',...
    'A_EW_spectra_max_period','A_EW_spectra_seconds02','A_EW_spectra_seconds05', ...
    'A_EW_spectra_seconds1', 'A_EW_spectra_seconds2', 'A_EW_spectra_seconds5',...
    'A_EW_HV_amp','A_EW_HV_freq',...
    'rotated_PGA', 'A_rotated_arias_intensity','A_rotated_arias_D595', ...
    'A_rotated_arias_D575', 'A_rotated_arias_rate', 'A_rotated_spectra_max',...
    'A_rotated_spectra_max_period','A_rotated_spectra_seconds02','A_rotated_spectra_seconds05', ...
    'A_rotated_spectra_seconds1', 'A_rotated_spectra_seconds2', 'A_rotated_spectra_seconds5',...
    'A_rotated_HV_amp','A_rotated_HV_freq', 'V_PGA', 'A_V_arias_intensity','A_V_arias_D595', ...
    'A_V_arias_D575', 'A_V_arias_rate', 'A_V_spectra_max',...
    'A_V_spectra_max_period','A_V_spectra_seconds02','A_V_spectra_seconds05', ...
    'A_V_spectra_seconds1', 'A_V_spectra_seconds2', 'A_V_spectra_seconds5',...
    'A_complex_HV_amp','A_complex_HV_freq',...
    'NS_PGV',  'V_NS_spectra_max', 'V_NS_spectra_max_period',...
    'V_NS_spectra_seconds_02', 'V_NS_spectra_seconds_05', 'V_NS_spectra_seconds_1'...,
    'V_NS_spectra_seconds_2', 'V_NS_spectra_seconds_5',...
    'EW_PGV',  'V_EW_spectra_max', 'V_EW_spectra_max_period',...
    'V_EW_spectra_seconds_02', 'V_EW_spectra_seconds_05', 'V_EW_spectra_seconds_1'...,
    'V_EW_spectra_seconds_2', 'V_EW_spectra_seconds_5',...
    'V_PGV',  'V_V_spectra_max', 'V_V_spectra_max_period',...
    'V_V_spectra_seconds_02', 'V_V_spectra_seconds_05', 'V_V_spectra_seconds_1'...,
    'V_V_spectra_seconds_2', 'V_V_spectra_seconds_5',...
    'rotated_PGV',  'V_rotated_spectra_max', 'V_rotated_spectra_max_period',...
    'V_rotated_spectra_seconds_02', 'V_rotated_spectra_seconds_05', 'V_rotated_spectra_seconds_1'...,
    'V_rotated_spectra_seconds_2', 'V_rotated_spectra_seconds_5',...
    'NS_PGD',  'D_NS_spectra_max', 'D_NS_spectra_max_period',...
    'D_NS_spectra_seconds_02', 'D_NS_spectra_seconds_05', 'D_NS_spectra_seconds_1', ...
    'D_NS_spectra_seconds_2', 'D_NS_spectra_seconds_5',...
    'EW_PGD',  'D_EW_spectra_max', 'D_EW_spectra_max_period',...
    'D_EW_spectra_seconds_02', 'D_EW_spectra_seconds_05', 'D_EW_spectra_seconds_1',...
    'D_EW_spectra_seconds_2', 'D_EW_spectra_seconds_5',...
    'V_PGD',  'D_V_spectra_max', 'D_V_spectra_max_period',...
    'D_V_spectra_seconds_02', 'D_V_spectra_seconds_05', 'D_V_spectra_seconds_1',...
    'D_V_spectra_seconds_2', 'D_V_spectra_seconds_5',...
    'rotated_PGD',  'D_rotated_spectra_max', 'D_rotated_spectra_max_period',...
    'D_rotated_spectra_seconds_02', 'D_rotated_spectra_seconds_05', 'D_rotated_spectra_seconds_1',...
    'D_rotated_spectra_seconds_2', 'D_rotated_spectra_seconds_5'};
for i = 2% : length(stationlist)
    station = stationlist(i).name;
    cd(strcat(datapath, '\', station))
    eventlist = dir;
    eventlist = eventlist(3:length(eventlist));
%     event_num{i,1} = station;
%     event_num{i,2} = length(eventlist);
    for j = 1: length(eventlist)
        counter = counter +1;
        disp(counter)
        filename = eventlist(j);
        filename = strcat(filename.folder, '\', filename.name);
        load(filename)
        % event
        A{counter, 1} = data.meta.event.ID; % event ID
        A{counter, 2} = data.meta.event.mag; % event magnitude
        A{counter, 3} = data.meta.event.lat; % event latitude
        A{counter, 4} = data.meta.event.lon; % event longitude
        A{counter, 5} = data.meta.event.depth; % event depth
        A{counter, 6} = data.meta.event.azimuth; % event azimuth
        A{counter, 7} = data.meta.event.epi_dist; % event epicentral distance
        
        % station
        A{counter, 8} = data.meta.station.name; % station name
        A{counter, 9} = data.meta.station.lat; % station latitude
        A{counter, 10} = data.meta.station.lon; % station longitude
        A{counter, 11} = data.meta.station.alt; % station altitude
        A{counter, 12} = data.meta.station.soil; % station soil
        
        % instrument
        A{counter, 13} = data.meta.instrument.model; % instrument model
        A{counter, 14} = data.meta.instrument.fs; % instrument frequency
        
        % ACCELERATION
        
        % NS
        % arias
        A{counter, 15} = data.processing.filtereddata.acceleration.NS.PGA; % PGA
        A{counter, 16} = data.processing.filtereddata.acceleration.NS.arias.intensity; % arias intensity
        A{counter, 17} = data.processing.filtereddata.acceleration.NS.arias.D595; % arias D595
        A{counter, 18} = data.processing.filtereddata.acceleration.NS.arias.D575; % arias D575
        A{counter, 19} = data.processing.filtereddata.acceleration.NS.arias.rate; % arias rate
        
        % response spectra
        A{counter, 20} = data.processing.filtereddata.acceleration.NS.spectra.max; % maximum spectra
        A{counter, 21} = data.processing.filtereddata.acceleration.NS.spectra.max_period; % period of maximum spectra
        A{counter, 22} = data.processing.filtereddata.acceleration.NS.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 23} = data.processing.filtereddata.acceleration.NS.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 24} = data.processing.filtereddata.acceleration.NS.spectra.seconds_1; % response at 1 second
        A{counter, 25} = data.processing.filtereddata.acceleration.NS.spectra.seconds_2; % response at 2 seconds
        A{counter, 26} = data.processing.filtereddata.acceleration.NS.spectra.seconds_5; % response at 5 seconds

        % HVSR
        A{counter, 27} = data.processing.filtereddata.acceleration.NS.HVSR.smooth.Amp; % amplification
        A{counter, 28} = data.processing.filtereddata.acceleration.NS.HVSR.smooth.fn; % frequency
        
        % EW
        % arias
        A{counter, 29} = data.processing.filtereddata.acceleration.EW.PGA; % PGA
        A{counter, 30} = data.processing.filtereddata.acceleration.EW.arias.intensity; % arias intensity
        A{counter, 31} = data.processing.filtereddata.acceleration.EW.arias.D595; % arias D595
        A{counter, 32} = data.processing.filtereddata.acceleration.EW.arias.D575; % arias D575
        A{counter, 33} = data.processing.filtereddata.acceleration.EW.arias.rate; % arias rate
        
        % response spectra
        A{counter, 34} = data.processing.filtereddata.acceleration.EW.spectra.max; % maximum spectra
        A{counter, 35} = data.processing.filtereddata.acceleration.EW.spectra.max_period; % period of maximum spectra
        A{counter, 36} = data.processing.filtereddata.acceleration.EW.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 37} = data.processing.filtereddata.acceleration.EW.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 38} = data.processing.filtereddata.acceleration.EW.spectra.seconds_1; % response at 1 second
        A{counter, 39} = data.processing.filtereddata.acceleration.EW.spectra.seconds_2; % response at 2 seconds
        A{counter, 40} = data.processing.filtereddata.acceleration.EW.spectra.seconds_5; % response at 5 seconds

        % HVSR
        A{counter, 41} = data.processing.filtereddata.acceleration.EW.HVSR.smooth.Amp; % amplification
        A{counter, 42} = data.processing.filtereddata.acceleration.EW.HVSR.smooth.fn; % frequency
        
        % rotated
        % arias
        A{counter, 43} = data.processing.filtereddata.acceleration.rotated.PGA; % PGA
        A{counter, 44} = data.processing.filtereddata.acceleration.rotated.arias.intensity; % arias intensity
        A{counter, 45} = data.processing.filtereddata.acceleration.rotated.arias.D595; % arias D595
        A{counter, 46} = data.processing.filtereddata.acceleration.rotated.arias.D575; % arias D575
        A{counter, 47} = data.processing.filtereddata.acceleration.rotated.arias.rate; % arias rate
        
        % response spectra
        A{counter, 48} = data.processing.filtereddata.acceleration.rotated.spectra.max; % maximum spectra
        A{counter, 49} = data.processing.filtereddata.acceleration.rotated.spectra.max_period; % period of maximum spectra
        A{counter, 50} = data.processing.filtereddata.acceleration.rotated.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 51} = data.processing.filtereddata.acceleration.rotated.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 52} = data.processing.filtereddata.acceleration.rotated.spectra.seconds_1; % response at 1 second
        A{counter, 53} = data.processing.filtereddata.acceleration.rotated.spectra.seconds_2; % response at 2 seconds
        A{counter, 54} = data.processing.filtereddata.acceleration.rotated.spectra.seconds_5; % response at 5 seconds

        % HVSR
        A{counter, 55} = data.processing.filtereddata.acceleration.rotated.HVSR.smooth.Amp; % amplification
        A{counter, 56} = data.processing.filtereddata.acceleration.rotated.HVSR.smooth.fn; % frequency
        
        % V
        % arias
        A{counter, 57} = data.processing.filtereddata.acceleration.V.PGA; % PGA
        A{counter, 58} = data.processing.filtereddata.acceleration.V.arias.intensity; % arias intensity
        A{counter, 59} = data.processing.filtereddata.acceleration.V.arias.D595; % arias D595
        A{counter, 60} = data.processing.filtereddata.acceleration.V.arias.D575; % arias D575
        A{counter, 61} = data.processing.filtereddata.acceleration.V.arias.rate; % arias rate
        
        % response spectra
        A{counter, 62} = data.processing.filtereddata.acceleration.V.spectra.max; % maximum spectra
        A{counter, 63} = data.processing.filtereddata.acceleration.V.spectra.max_period; % period of maximum spectra
        A{counter, 64} = data.processing.filtereddata.acceleration.V.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 65} = data.processing.filtereddata.acceleration.V.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 66} = data.processing.filtereddata.acceleration.V.spectra.seconds_1; % response at 1 second
        A{counter, 67} = data.processing.filtereddata.acceleration.V.spectra.seconds_2; % response at 2 seconds
        A{counter, 68} = data.processing.filtereddata.acceleration.V.spectra.seconds_5; % response at 5 seconds
        
        % complex
        
        % HVSR
        A{counter, 69} = data.processing.filtereddata.acceleration.complex.HVSR.smooth.Amp; % amplification
        A{counter, 70} = data.processing.filtereddata.acceleration.complex.HVSR.smooth.fn; % frequency
        
        % VELOCITY
        
        % NS
        A{counter, 71} = data.processing.filtereddata.velocity.NS.PGV; % PGV
        
        % spectra
        A{counter, 72} = data.processing.filtereddata.velocity.NS.spectra.max; % maximum velocity spectra
        A{counter, 73} = data.processing.filtereddata.velocity.NS.spectra.max_period; % period of maximum spectra
        A{counter, 74} = data.processing.filtereddata.velocity.NS.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 75} = data.processing.filtereddata.velocity.NS.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 76} = data.processing.filtereddata.velocity.NS.spectra.seconds_1; % response at 1 second
        A{counter, 77} = data.processing.filtereddata.velocity.NS.spectra.seconds_2; % response at 2 seconds
        A{counter, 78} = data.processing.filtereddata.velocity.NS.spectra.seconds_5; % response at 5 seconds

        % EW
        A{counter, 79} = data.processing.filtereddata.velocity.EW.PGV; % PGV
        
        % spectra
        A{counter, 80} = data.processing.filtereddata.velocity.EW.spectra.max; % maximum velocity spectra
        A{counter, 81} = data.processing.filtereddata.velocity.EW.spectra.max_period; % period of maximum spectra
        A{counter, 82} = data.processing.filtereddata.velocity.EW.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 83} = data.processing.filtereddata.velocity.EW.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 84} = data.processing.filtereddata.velocity.EW.spectra.seconds_1; % response at 1 second
        A{counter, 85} = data.processing.filtereddata.velocity.EW.spectra.seconds_2; % response at 2 seconds
        A{counter, 86} = data.processing.filtereddata.velocity.EW.spectra.seconds_5; % response at 5 seconds
        
        % V
        A{counter, 87} = data.processing.filtereddata.velocity.V.PGV; % PGV
        
        % spectra
        A{counter, 88} = data.processing.filtereddata.velocity.V.spectra.max; % maximum velocity spectra
        A{counter, 89} = data.processing.filtereddata.velocity.V.spectra.max_period; % period of maximum spectra
        A{counter, 90} = data.processing.filtereddata.velocity.V.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 91} = data.processing.filtereddata.velocity.V.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 92} = data.processing.filtereddata.velocity.V.spectra.seconds_1; % response at 1 second
        A{counter, 93} = data.processing.filtereddata.velocity.V.spectra.seconds_2; % response at 2 seconds
        A{counter, 94} = data.processing.filtereddata.velocity.V.spectra.seconds_5; % response at 5 seconds
       
        % rotated
        A{counter, 95} = data.processing.filtereddata.velocity.rotated.PGV; % PGV
        
        % spectra
        A{counter, 96} = data.processing.filtereddata.velocity.rotated.spectra.max; % maximum velocity spectra
        A{counter, 97} = data.processing.filtereddata.velocity.rotated.spectra.max_period; % period of maximum spectra
        A{counter, 98} = data.processing.filtereddata.velocity.rotated.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 99} = data.processing.filtereddata.velocity.rotated.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 100} = data.processing.filtereddata.velocity.rotated.spectra.seconds_1; % response at 1 second
        A{counter, 101} = data.processing.filtereddata.velocity.rotated.spectra.seconds_2; % response at 2 seconds
        A{counter, 102} = data.processing.filtereddata.velocity.rotated.spectra.seconds_5; % response at 5 seconds
       
        
        % DISPLACEMENT
        
        % NS
        A{counter, 103} = data.processing.filtereddata.displacement.NS.PGD; % PGD    
        
        % spectra
        A{counter, 104} = data.processing.filtereddata.displacement.NS.spectra.max; % maximum displacement spectra
        A{counter, 105} = data.processing.filtereddata.displacement.NS.spectra.max_period; % period of maximum spectra
        A{counter, 106} = data.processing.filtereddata.displacement.NS.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 107} = data.processing.filtereddata.displacement.NS.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 108} = data.processing.filtereddata.displacement.NS.spectra.seconds_1; % response at 1 second
        A{counter, 109} = data.processing.filtereddata.displacement.NS.spectra.seconds_2; % response at 2 seconds
        A{counter, 110} = data.processing.filtereddata.displacement.NS.spectra.seconds_5; % response at 5 seconds
    
        % EW
        A{counter, 111} = data.processing.filtereddata.displacement.EW.PGD; % PGD    
        
        % spectra
        A{counter, 112} = data.processing.filtereddata.displacement.EW.spectra.max; % maximum displacement spectra
        A{counter, 113} = data.processing.filtereddata.displacement.EW.spectra.max_period; % period of maximum spectra
        A{counter, 114} = data.processing.filtereddata.displacement.EW.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 115} = data.processing.filtereddata.displacement.EW.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 116} = data.processing.filtereddata.displacement.EW.spectra.seconds_1; % response at 1 second
        A{counter, 117} = data.processing.filtereddata.displacement.EW.spectra.seconds_2; % response at 2 seconds
        A{counter, 118} = data.processing.filtereddata.displacement.EW.spectra.seconds_5; % response at 5 seconds    
    
        % V
        A{counter, 119} = data.processing.filtereddata.displacement.V.PGD; % PGD    
        
        % spectra
        A{counter, 120} = data.processing.filtereddata.displacement.V.spectra.max; % maximum displacement spectra
        A{counter, 121} = data.processing.filtereddata.displacement.V.spectra.max_period; % period of maximum spectra
        A{counter, 122} = data.processing.filtereddata.displacement.V.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 123} = data.processing.filtereddata.displacement.V.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 124} = data.processing.filtereddata.displacement.V.spectra.seconds_1; % response at 1 second
        A{counter, 125} = data.processing.filtereddata.displacement.V.spectra.seconds_2; % response at 2 seconds
        A{counter, 126} = data.processing.filtereddata.displacement.V.spectra.seconds_5; % response at 5 seconds    
    
        % rotated
        A{counter, 127} = data.processing.filtereddata.displacement.rotated.PGD; % PGD    
        
        % spectra
        A{counter, 128} = data.processing.filtereddata.displacement.rotated.spectra.max; % maximum displacement spectra
        A{counter, 129} = data.processing.filtereddata.displacement.rotated.spectra.max_period; % period of maximum spectra
        A{counter, 130} = data.processing.filtereddata.displacement.rotated.spectra.seconds_02; % response at 0.2 seconds
        A{counter, 131} = data.processing.filtereddata.displacement.rotated.spectra.seconds_05; % response at 0.5 seconds
        A{counter, 132} = data.processing.filtereddata.displacement.rotated.spectra.seconds_1; % response at 1 second
        A{counter, 133} = data.processing.filtereddata.displacement.rotated.spectra.seconds_2; % response at 2 seconds
        A{counter, 134} = data.processing.filtereddata.displacement.rotated.spectra.seconds_5; % response at 5 seconds    
    
    end
end