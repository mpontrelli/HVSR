% Program to convert one day of Amesbury Town Hall Data from Reftek to SAC
% data format

% This program creates a batch file in text format that can be copied and
% pasted into a terminal command window.  The program will work from any
% folder in the terminal command window.  When executed, it should produce
% a sac subdirectory in the 1 folder where the sac files are stored.
% Sometimes it skips reading one or more of the data files to be converted
% to sac.  If this happens, paste the batch file contents into the command
% window again, and it should find the remaining files and convert them to
% sac format (if not, paste it in again.

clear all

% Begin by specifying the year, and day of the data to be converted
%Yearstr='2009';
%Daystr='007';
Yearstr=input('Give year: ','s');
Daystr=input('Give day: ','s');

% % %this allows you to run function as if it is a dos command%%%
% [status,cmdout] = dos( 'pas2sac filename' );
% status, cmdout
% [status, cmdout] = dos('dir');
% status, cmdout

%%%combine all 3 into 1
% Waterboro data section
thispath{1}='/Users/Alex/Downloads/converters/macos64/pas2sac/';
thisdir{1}=[thispath{1},Yearstr,Daystr,'/9B27/1/000000004_0036EE88/'];
thispath{2}='/Users/ebel/Documents/WaterboroME_EQ/Sta2McCardleyData/';
thisdir{2}=[thispath{2},Yearstr,Daystr,'/B680/1'];
thispath{3}='/Users/ebel/Documents/WaterboroME_EQ/Sta3KellyData/';
thisdir{3}=[thispath{3},Yearstr,Daystr,'/9B27/1'];

%%%This needs to open PAS2SAC instead%%%
fiddir2=fopen('runreftek.bat','w');
tstfl=fclose(fiddir2);


%%Write batch commands into here and run that-- OR use a dos command in
%%the future

for jjtt=1:3
fiddir2=fopen('runreftek.bat','a');
fprintf(fiddir2,'cd %s\n',thisdir{jjtt});
tstfl=fclose(fiddir2);


% Here we get a listing of files from the subdirectory with the data
fiddir=fopen('datadir.bat','w');
fprintf(fiddir,'#!/bin/bash\n');
fprintf(fiddir,'ls -a %s > refdatadir.lst\n',thisdir{jjtt});
tstfl=fclose(fiddir);

% Next execute this batch file
if tstfl == 0
    !chmod a=rwx datadir.bat
    !./datadir.bat
end  % if tstfl == 0



% We next open the output file that was just created and read the list
% of station data to be processed
%ON A UNIX%, refdatadir is output
fidst=fopen('refdatadir.lst');
jsta=0;
if fidst>0
    % Read the name of the next file and process it if it is a Reftek data
    %  file
%     Again, can probably consolidate sections of this %%%
    fstatest='1';
    while fstatest~=' '
        fstaline=fgetl(fidst);  % Read the next line in the file
        
        % The next couple of if statements check to see if we have the name
        % of a real data file to be processed
        if fstaline == -1  % If this is true, then we have reached the end
               % of the data file
            fstatest=' ';
        end  % if fstaline == -1
        linelen=length(fstaline);
        if linelen > 10 && strcmp(fstaline(10:10),'_')==1
            nextsta=fstaline; % If we were true to this point, then we
              % have a real data file to be processed
              %%%Should we keep this as linelen>10???
        else
            nextsta='';
        end  % if linelen > 10 && strcmp(fstaline(10:10),'_')==1
        
        % Next, if we have a real data file, then process it
        if ~isempty(nextsta)
            thisfile2=[thisdir{jjtt},'/',fstaline];
            fiddir2=fopen('runreftek.bat','a');
            fprintf(fiddir2,'ref2segy -f %s\n',thisfile2);
            %fprintf(fiddir2,'100\n');  % For the Amesbury data
            fprintf(fiddir2,'40\n');  % For the Waterboro data
            fprintf(fiddir2,'32\n');
            fprintf(fiddir2,'1-3\n');
            fprintf(fiddir2,'1\n');
            fprintf(fiddir2,'1\n');
            fprintf(fiddir2,'1\n');
            fprintf(fiddir2,'\n\n\n\n\n\n\n\n\n\n');
     
            %%%just delete all of this???
            %fprintf(fiddir2,'ref2segy -f %s\n',thisfile2);
            tstfl=fclose(fiddir2);
            % Next execute this batch file
            if tstfl == 0
               !chmod a=rwx runreftek.bat
               !./runreftek.bat
            end  % if tstfl == 0
            %%% Maybe above piece is NOT needed???
            
            %fstatest=' ';
            
        end  % if length(nextsta)>0
        
        
            
    end  % while fstatest~=' '
end  % if fidst>0

% We now have all of the data translated to segy file format.  Here we go
% into the subdirectory and convert those files to SAC format
Rfile=['R',Daystr,'.01'];
fiddir2=fopen('runreftek.bat','a');
fprintf(fiddir2,'cd %s\n',Rfile);
fprintf(fiddir2,'/Users/ebel/Documents/Matlab/segy2sac *.*\n');
tstfl=fclose(fiddir2);

fclose(fidst);

end  % for jjtt=1:3