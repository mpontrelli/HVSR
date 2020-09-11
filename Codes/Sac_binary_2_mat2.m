%% Read sac files and turn them into matlab files
Vfname = 'C:\Users\Marshall\Box\Geohazards Research Group\2020_7_1\Charles_test\3c\2020183180806005_pembr_1_1.sac';
NSfname = 'C:\Users\Marshall\Box\Geohazards Research Group\2020_7_1\Charles_test\3c\2020183180806005_pembr_1_2.sac';
EWfname = 'C:\Users\Marshall\Box\Geohazards Research Group\2020_7_1\Charles_test\3c\2020183180806005_pembr_1_3.sac';
outpath = 'C:\Users\Marshall\Box\Geohazards Research Group\2020_7_1\Charles_test\3c';
name = '3c';

%% first file
[V] = ReadSacBinaryFile(Vfname); %vertical
[NS] = ReadSacBinaryFile(NSfname); %North-south
[EW] = ReadSacBinaryFile(EWfname); %East-West

%% Second file if you need it 
Vfname2 = 'C:\Users\Marshall\Box\Geohazards Research Group\2020_7_1\Charles_test\3b\2020183180806005_pembr_1_1.sac';
NSfname2 = 'C:\Users\Marshall\Box\Geohazards Research Group\2020_7_1\Charles_test\3b\2020183180806005_pembr_1_2.sac';
EWfname2 = 'C:\Users\Marshall\Box\Geohazards Research Group\2020_7_1\Charles_test\3b\2020183180806005_pembr_1_3.sac';

[V2] = ReadSacBinaryFile(Vfname2); %vertical
[NS2] = ReadSacBinaryFile(NSfname2); %North-south
[EW2] = ReadSacBinaryFile(EWfname2); %East-West


%% Now concatenate
V = vertcat(V,V2);
NS = vertcat(NS,NS2);
EW = vertcat(EW, EW2);

data.V = V;
data.NS = NS;
data.EW = EW;

save(strcat(outpath, '\',name, '.mat'), 'data')