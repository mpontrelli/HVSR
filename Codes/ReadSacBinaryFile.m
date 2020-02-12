function [waveform, delta, stla, stlo] = ReadSacBinaryFile(filename)
% ReadSacBinary.m
%
% Matlab m-file to read Binary SAC files
%
% This file is currently set up to read only evenly spaced waveform files.
%
% First, get the name of the file
filnam=filename;
%filnam=datafilename;

%  Open the file
% [fid,message]=fopen(filnam,'r');
[fid,message]=fopen(filnam,'r','b');  % The 'b' argument is for a big-endian file
% Note:  Becaus e SAC files were designed for SUN computers, which are
%  "big-endian" machines, all SAC files are written with the binary data
%  in big-endian format.  On a PC, the natural way that numbers are stored
%  is in little-endian format.  Thus, when a SAC file is opened on a PC,
%  the "open" statement must include a qualifier indicating that the binary
%  data file was written in big-endian format.  If this macro is being run
%  on a big-endian machine (like a Sun, Mac, HP, IBM RS or SGI), then the
%  'b' argument in the open statement is optional.
%
if fid == -1  % Oops, there was an error and no file was opened
	info='Cannot open file:';
	filnam;
    message;
  else
    % We successfully opened the file, so read it
    
    % Parse out all the SAC header variables
    [delta,count]=fread(fid,1,'float32');
    [depmin,count]=fread(fid,1,'float32');
    [depmax,count]=fread(fid,1,'float32');
    [scale,count]=fread(fid,1,'float32');
    [odelta,count]=fread(fid,1,'float32');
    [b,count]=fread(fid,1,'float32');
    [e,count]=fread(fid,1,'float32');
    [o,count]=fread(fid,1,'float32');
    [a,count]=fread(fid,1,'float32');
    [junk1,count]=fread(fid,1,'float32');
    [t0,count]=fread(fid,1,'float32');
    [t1,count]=fread(fid,1,'float32');
    [t2,count]=fread(fid,1,'float32');
    [t3,count]=fread(fid,1,'float32');
    [t4,count]=fread(fid,1,'float32');
    [t5,count]=fread(fid,1,'float32');
    [t6,count]=fread(fid,1,'float32');
    [t7,count]=fread(fid,1,'float32');
    [t8,count]=fread(fid,1,'float32');
    [t9,count]=fread(fid,1,'float32');
    [f,count]=fread(fid,1,'float32');
    [resp0,count]=fread(fid,1,'float32');
    [resp1,count]=fread(fid,1,'float32');
    [resp2,count]=fread(fid,1,'float32');
    [resp3,count]=fread(fid,1,'float32');
    [resp4,count]=fread(fid,1,'float32');
    [resp5,count]=fread(fid,1,'float32');
    [resp6,count]=fread(fid,1,'float32');
    [resp7,count]=fread(fid,1,'float32');
    [resp8,count]=fread(fid,1,'float32');
    [resp9,count]=fread(fid,1,'float32');
    [stla,count]=fread(fid,1,'float32');
    [stlo,count]=fread(fid,1,'float32');
    [stel,count]=fread(fid,1,'float32');
    [stdp,count]=fread(fid,1,'float32');
    [evla,count]=fread(fid,1,'float32');
    [evlo,count]=fread(fid,1,'float32');
    [evel,count]=fread(fid,1,'float32');
    [evdp,count]=fread(fid,1,'float32');
    [junk2,count]=fread(fid,1,'float32');
    [user0,count]=fread(fid,1,'float32');
    [user1,count]=fread(fid,1,'float32');
    [user2,count]=fread(fid,1,'float32');
    [user3,count]=fread(fid,1,'float32');
    [user4,count]=fread(fid,1,'float32');
    [user5,count]=fread(fid,1,'float32');
    [user6,count]=fread(fid,1,'float32');
    [user7,count]=fread(fid,1,'float32');
    [user8,count]=fread(fid,1,'float32');
    [user9,count]=fread(fid,1,'float32');
    [dist,count]=fread(fid,1,'float32');
    [az,count]=fread(fid,1,'float32');
    [baz,count]=fread(fid,1,'float32');
    [gcarc,count]=fread(fid,1,'float32');
    [junk3,count]=fread(fid,1,'float32');
    [junk4,count]=fread(fid,1,'float32');
    [depmen,count]=fread(fid,1,'float32');
    [cmpaz,count]=fread(fid,1,'float32');
    [cmpinc,count]=fread(fid,1,'float32');
    [junk5,count]=fread(fid,1,'float32');
    [junk6,count]=fread(fid,1,'float32');
    [junk7,count]=fread(fid,1,'float32');
    [junk8,count]=fread(fid,1,'float32');
    [junk9,count]=fread(fid,1,'float32');
    [junk10,count]=fread(fid,1,'float32');
    [junk11,count]=fread(fid,1,'float32');
    [junk12,count]=fread(fid,1,'float32');
    [junk13,count]=fread(fid,1,'float32');
    [junk14,count]=fread(fid,1,'float32');
    [junk15,count]=fread(fid,1,'float32');
    [nzyear,count]=fread(fid,1,'uint32');
    [nzjday,count]=fread(fid,1,'uint32');
    [nzhour,count]=fread(fid,1,'uint32');
    [nzmin,count]=fread(fid,1,'uint32');
    [nzsec,count]=fread(fid,1,'uint32');
    [nzmsec,count]=fread(fid,1,'uint32');
    [nvhdr,count]=fread(fid,1,'uint32');
    [junk16,count]=fread(fid,1,'uint32');
    [junk17,count]=fread(fid,1,'uint32');
    [npts,count]=fread(fid,1,'uint32');
    [junk100,count]=fread(fid,1,'uint32');
    [junk18,count]=fread(fid,1,'uint32');
    [junk19,count]=fread(fid,1,'uint32');
    [junk20,count]=fread(fid,1,'uint32');
    [junk21,count]=fread(fid,1,'uint32');
    [iftype,count]=fread(fid,1,'uint32');
    [idep,count]=fread(fid,1,'uint32');
    [iztype,count]=fread(fid,1,'uint32');
    [junk22,count]=fread(fid,1,'uint32');
    [iinst,count]=fread(fid,1,'uint32');
    [istreg,count]=fread(fid,1,'uint32');
    [ievreg,count]=fread(fid,1,'uint32');
    [ievtyp,count]=fread(fid,1,'uint32');
    [iqual,count]=fread(fid,1,'uint32');
    [isynth,count]=fread(fid,1,'uint32');
    [junk23,count]=fread(fid,1,'uint32');
    [junk24,count]=fread(fid,1,'uint32');
    [junk25,count]=fread(fid,1,'uint32');
    [junk26,count]=fread(fid,1,'uint32');
    [junk27,count]=fread(fid,1,'uint32');
    [junk28,count]=fread(fid,1,'uint32');
    [junk29,count]=fread(fid,1,'uint32');
    [junk30,count]=fread(fid,1,'uint32');
    [junk31,count]=fread(fid,1,'uint32');
    [junk32,count]=fread(fid,1,'uint32');
    [leven,count]=fread(fid,1,'uint32');
    [lspol,count]=fread(fid,1,'uint32');
    [lovrok,count]=fread(fid,1,'uint32');
    [lcalda,count]=fread(fid,1,'uint32');
    [junk33,count]=fread(fid,1,'uint32');
    [ktemp,count]=fread(fid,8,'schar');
    kstnm=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kevnm=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kjunk=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    khole=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    ko=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    ka=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kt0=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kt1=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kt2=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kt3=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kt4=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kt5=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kt6=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kt7=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kt8=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kt9=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kf=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kuser0=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kuser1=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kuser2=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kcmpnm=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    knetwk=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kdatrd=char(ktemp');
    [ktemp,count]=fread(fid,8,'schar');
    kinst=char(ktemp');
    % Done reading the header, now read the data
    [waveform,countw]=fread(fid,inf,'float32');
end  %  End the if loop
fclose(fid);  % Close the file
%
% Send some basic info to the Matlab Command window
info = ['Station = ',kstnm,',  Comp = ',kcmpnm];
fprintf('%s\n',info);
info = ['Number of points = ',num2str(npts),',  DT = ',num2str(delta)];
fprintf('%s\n',info);
info = 'The seismogram is in the array: "waveform"';