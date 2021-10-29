%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        HV-Inv 1.31 Beta GUI
%
% Inversion of the H/V spectral ratio based on the diffuse wavefield theory
%
% Matlab functions deal with GUIs and inversion algorithms. HVTI is the
% main program of HV-Inv 1.3 Beta.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Copyright (C) 2014,2016 José Piña-Flores, Antonio García-Jerez.
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License version 3 as
%     published by the Free Software Foundation.
%
%     This program is distributed in the hope that it will be useful, but
%     WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.

function varargout = HVTI(varargin)
%% Function for HV-Inv GUI
% Type of OS
if ispc % Windows
    addpath(genpath('./bin/Functions/Func_Fig/Fig_Win'));
elseif ismac %Mac
    addpath(genpath('./bin/Functions/Func_Fig/Fig_Mac'));
elseif isunix %Linux
    addpath(genpath('./bin/Functions/Func_Fig/Fig_Linux'));
end
gui_Singleton = 1;
gui_State = struct('gui_Name',mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @HVTI_OpeningFcn, ...
    'gui_OutputFcn',  @HVTI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
%#ok<*DEFNU>
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      Variables for GUI and  BASE workspace      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Assign default values starting the GUI
function HVTI_OpeningFcn(hObject, ~, handles,  ~)
img = imread('LogoHVinv.png');
%% Display logo Java
jimg = im2java(img);
frame = javax.swing.JFrame;
frame.setUndecorated(true);
icon = javax.swing.ImageIcon(jimg);
label = javax.swing.JLabel(icon);
frame.getContentPane.add(label);
frame.pack;
imgSize = size(img);
frame.setSize(imgSize(2),imgSize(1));
screenSize = get(0,'ScreenSize');  % Get the screen size from the root object
frame.setLocation((screenSize(3)-imgSize(2))/2,...  % Center on the screen
    (screenSize(4)-imgSize(1))/2);
frame.show;
pause(1)
frame.hide;
if ismac
    setenv ('DYLD_LIBRARY_PATH','/usr/local/bin:/opt/local/lib:')
end
% GUI Variables
handles.output = hObject;
guidata(hObject, handles);
%% Kernels control
FGP=parallel(1);
if FGP==1
    handles.PARALL.Enable= 'off';
    open=2;% parallel computing not available
else
    LoadV('NUMW',FGP);
    open=0;% Parallel computing is currently dissabled. (=1 means enabled)
end
%% Variables for inversion
% Iterations,% Number of iterations % Number of last iterations,Search radius (for MCS, SA and MSA)
[A.Nsubiter,A.Niter,A.Niterf,A.NN]=deal(0,100,100,5); % Number of temperatures
%% Assignment of variables to the base workspace
%% Hard half space: Must the halfspace have the highest Vs?
%HHS  0 = no; 1 = yes HSS=1 default
%% Initial model
% ini1=[]
%% Identificator  of graphical panels in GUI
%'Perr',[] Global graphic error
%'fgd1',[] Legend of graphic of error
%'F1',[] Graphic of misfit
% 'F2',[]Graphic of current misfit
%INPO=50 Default for the number of models in the random initial population:
%% Tinv is the code for the inversion algorithm. MCS is the default method
%1=Monte Carlo Sampling (MCS), 2=Simulated Annealing (SA), 3=Modified SA (MSA), 4=Simplex ,5=Interior point (IP)
%'Tinv=1 Monte Carlo Samplig
% Number of layers NC=0;
%% Default settings for calling the forward H/V program:
% Surface waves
% 'ramasR'=5 Max. number of Rayleigh modes
% 'ramasL',5 Max. number of Love modes
% Body wave integrals
%'kb_by_dk',1000 Minimum points integration
%'mkb_by_dk',2000 Maximum points integration
% Low S-wave velocity zones (Vs decreases as depth increases for some
% layers) 'ZLVs'=0);% 0 = forbidden; 1 = allowed
% 0 = forbidden; 1 = allowed
% Low P-wave velocity zones (Vs decreases as depth increases for some
% layers)'ZLVp'=1); 0 = forbidden; 1 = allowed
% 0 = forbidden; 1 = allowed
% Initial Temperature for SA and MSA
% 1 = use a method based on a probability of increasing the misfit of the initial model.
% If T0=1, the probability of acceptance (in the begining) of a model which
% multiplies the misfit L2 of the initial model by EA is AP, that is,
% exp(-(EA*L2)/T)=AP or T=-A.EA*A.L2/log(A.AP)
%'T0'=1
%% Temperature
%'EA'=0.1= Relative misfit increment
%  'AP',0.5 (50%)Acceptance probability
%'RT',0.9 Reduction factor for Temperature
%% Stabilization factor (slight damping for body wave integrals)
% apsv'=0.005 ;0 = no damping

%% Load variables at workspace
LoadV('Tinv',1,'F2',[],'F1',[],'fgd1',[],'Perr',[],'ini1',[],'A',A,'open',open,'HHS',1,'INPO',50,'NC',0,'kb_by_dk',1000,'mkb_by_dk',2000,...
    'ramasL',5,'ramasR',5,'ZLVs',0,'ZLVp',1,'T0',1,'EA',0.1,'RT',0.9,'AP',0.5,'apsv',0.005);
%% Graphic control of the GUI
cla(handles.ploterror,'reset')
[handles.START.Visible,handles.edit1.Enable,handles.PS.Enable,handles.Models.Enable,...
    handles.Pon.Checked,handles.LP.Enable,handles.SPara.Enable,handles.WPara.Enable,...
    handles.ploterror.Visible,handles.LZVS.Checked,handles.MSA.Checked,handles.SA.Checked,...
    handles.MW.Checked,handles.LS.Checked,handles.LF.Checked,handles.LV.Enable,handles.Inversion.Enable,...
    handles.T0.Enable,handles.Redu.Enable,handles.Inversiont.Enable]=deal('off');
[handles.HHS.Checked,handles.T0misfit.Checked,handles.Poff.Checked,handles.LZVP.Checked]=deal('on');
%     [handles.T0.Enable,handles.Redu.Enable,handles.Inversiont.Enable]=deal('off');

%% GUI Menu Control
Open = findall(handles.figure1, 'tag', 'Standard.FileOpen');
delete(Open)
Open = findall(handles.figure1, 'tag', 'Exploration.Rotate');
delete(Open)
Open = findall(handles.figure1, 'tag', 'Exploration.Brushing');
delete(Open)
Open = findall(handles.figure1, 'tag',  'DataManager.Linking');
delete(Open)
Open = findall(handles.figure1, 'tag',  'Annotation.InsertLegend');
delete(Open)
Open = findall(handles.figure1, 'tag',  'Annotation.InsertColorbar');
delete(Open)
Open = findall(handles.figure1, 'tag',  'Standard.SaveFigure');
delete(Open)
Open = findall(handles.figure1, 'tag', 'Standard.NewFigure');
delete(Open)
Open = findall(handles.figure1, 'tag', 'Standard.EditPlot');
delete(Open)
movegui(handles.figure1,'center')

%% Control of the table of model parameters ranges(GUI)

function edit1_Callback (hObject, ~, handles)
%% Set number of layers (halfspace included)
NC = str2double(hObject.String);
if isnan(NC) || isempty(NC) || isinf(NC)
    [hObject.String,NC]=deal( 2);
    errordlg('Input must be a number','Error');
elseif NC<=1 || not(mod(NC,1))==0
    [hObject.String,NC]=deal( 2);
    errordlg('The number of layers should be greater or equal to 2 and integer','Error');
end
[ini1,INPO]=GetV('ini1','INPO');
% Verifica si el modelo inicial tiene el mismo numero de capas
if NC==length(ini1) || INPO~=0;
    NCok=1;% Variable de verificacion del modelo inicial con la tabla de parametros
else
    choice = questdlg(sprintf('The number of layers does not match the initial model.\n\n Initial model will be deleted'),'Hola Antonio','Ok','Cancel','Ok');
    % Handle response
    switch choice
        case 'Ok'
            NCok=1;
        case 'Cancel'
            NCok=0;
            NC=GetV('NC');
            hObject.String=NC;
    end
end
%% Control of fields in the GUI (table of model parameters)
if NCok
    [handles.START.Visible,handles.TableParameters.Enable,...
        handles.SPara.Enable,handles.LV.Enable,handles.T0.Enable,...
        handles.Redu.Enable,handles.Inversiont.Enable]=deal('on');
    %% Initial ranges for the model parameters
    [data,var]=deal(zeros(NC,10),zeros(NC,15));% contains the table
    [data(:,1), data(:,2),data(:,3),data(:,4), data(:,5),data(:,6),...
        data(:,9), data(:,10)]  =deal(10,50,400,6000,200,3400,0.25,0.4);%default parameters H Vp Vs Density Poisson
    [data(:,7),data(:,8)]=deal(2000);% min and max density
    [data(end,1),data(end,2)]=deal(0);
    var(:,1:3:13)=(data(:,1:2:9));% Minima to matrix var
    var(:,2:3:14)=(data(:,2:2:10));% Maxima to matrix var
    var(:,3:3:15)=var(:,2:3:14)-var(:,1:3:13);
    limdd=sum(var(:,2)); %max. depth down to the halfspace
    %% Asignation of variables to workspace
    LoadV('NC',NC,'var',var,'data',data,'limdd',limdd,'INPO',50,'ini1',[])
    %Edicion GUI
    handles.TableParameters.Data=num2cell(data);
    % handles.TableParameters.ColumnEditable',true(1,10));
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% "Target" Structure %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Target is HV
function THV_Callback(~, ~, handles)
% Reading input file (HV only)
pass=LoadCurve(handles,1);
handles.WPara.Enable= 'on';
LoadV('flag',1)% flag=1 for inversion of H/V only
if pass
    nm=GetV('C.nmHV');% obteniendo numero de muestras de HV cargado
    if nm>100
        warndlg(sprintf(['The effective number of samples for HV curve is greater than 100 '...
            '(' num2str(nm) ').\n\n',...
            'This is likely to slow down the inversion process.\n\n',...
            'A usual number of samples is 50. Resample the curve if necessary.']),'!! Warning !!')
    end
    handles.Inversion.Enable='on';
%     [handles.T0.Enable,handles.Redu.Enable,handles.Inversiont.Enable]=deal('off');
end

%% Target DC
function TDC_Callback(~, ~, handles)
% Reading input file (Dispersion curve only)
pass=LoadCurve(handles,2);
handles.WPara.Enable='off';
LoadV('flag',2)% flag=2 for inversion of DC only
if pass
    nm=GetV('C.nmDC');
    if nm>100
        warndlg(sprintf(['The effective number of samples for dispersion curve is greater than 100 '...
            '(' num2str(nm) ').\n\n',...
            'This is likely to slow down the inversion process.\n\n',...
            'A usual number of samples is 50. Resample the curve if necessary.']),'!! Warning !!')
    end
    handles.Inversion.Enable='on';
%     [handles.T0.Enable,handles.Redu.Enable,handles.Inversiont.Enable]=deal('off');
end

%% Target HV+DC
function THVDC_Callback(~,~, handles)
% Reading input file (H/V and dispersion curve)
pass=LoadCurve(handles,3);
LoadV('flag',3) % flag=3 for joint inversion
handles.WPara.Enable='on';
if pass
    [nmD,nmH]=GetV('C.nmDC','C.nmHV');
    if nmH>100 && nmD>100
        warndlg(sprintf(['The effective number of samples for dispersion curve and DC curve ', '(' num2str(nmD) ') & HV curve ' '(' num2str(nmH) ')' ...
            ' are greater than 100.\n\n',...
            'This is likely to slow down the inversion process.\n\n',...
            'A usual number of samples is 50. Resample all curves if necessary.']),'!! Warning !!')
    elseif nmD>100
        warndlg(sprintf(['The effective number of samples for dispersion curve is greater than 100 '...
            '(' num2str(nmD) ').\n\n',...
            'This is likely to slow down the inversion process.\n\n',...
            'A usual number of samples is 50. Resample the curve if necessary.']),'!! Warning !!')
    elseif nmH>100
        warndlg(sprintf(['The effective number of samples for HV curve is greater than 100 '...
            '(' num2str(nmH) ').\n\n',...
            'This is likely to slow down the inversion process.\n\n',...
            'A usual number of samples is 50. Resample the curve if necessary.']),'!! Warning !!')
    end
    handles.Inversion.Enable='on';
%     [handles.T0.Enable,handles.Redu.Enable,handles.Inversiont.Enable]=deal('off');
    
end




%% Foward HV
% New GUI for forward problem
function FHV_Callback(~, ~, ~)
FHV;


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%Start Inversion%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function START_Callback(~, ~, handles)
[data,var,varoriginal,ZLVs,ZLVp,NC,HHS] =GetV ('data','var','var','ZLVs','ZLVp','NC','HHS');
[id,fow]=deal(0);%id error type in parameters table;Fow is id for update parameters
[cap2,cap1]=deal(cell(0));
nanok=sum(sum(isnan(var)));
if (sum(var(:,end-1)>=0.5)+sum(var(:,end-2)>=0.5))==0;
    if nanok==0
        [hok]=find(var(:,2)<varoriginal(:,1));
        if~isempty(hok)
            errordlg(sprintf('The minimum thickness values are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            id=3;
        end
        [vpok]=find(var(:,5)<varoriginal(:,4));
        if~isempty(vpok)
            errordlg(sprintf('The minimum vp values are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            id=3;
        end
        [vsok]=find(var(:,8)<varoriginal(:,7));
        if~isempty(vsok)
            errordlg(sprintf('The minimum vs values are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            id=3;
        end
        [rhook]=find(var(:,11)<varoriginal(:,10));
        if~isempty(rhook)
            errordlg(sprintf('The minimum density values are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            id=3;
        end
        [pok]=find(var(:,14)<varoriginal(:,13));
        if~isempty(pok)
            errordlg(sprintf('The minimum Poisson ratio values are larger than the maximum\n\n Please, change parameters'),'Bad Condition');
            id=3;
        end
        if id~=3
            if ZLVp~=0
                var(:,7)=max(var(:,7),var(:,4).*sqrt((1-(2*var(:,14)))./(2*(1-var(:,14)))));% We increase minimum Vs if is smaller than the result from maximum poi and minimum Vp
                var(:,8)=min(var(:,8),var(:,5).*sqrt((1-(2*var(:,13)))./(2*(1-var(:,13)))));% We decrease maximum Vs if is greater than the result from minimum poi and maximum Vp
                var(:,9)=var(:,8)-var(:,7);
            end
            if ZLVs==0
                % Replacing var(kl,8) (max Vs) with the minimum among the
                % maximum value of Vs in the kl-th layer and the maximun
                % values of the shalower layers
                for kl=NC-1:-1:1
                    var(kl,8)=min(var(kl,8),var(kl+1,8));
                end
                % Replacing var(kl,7) (min Vs) with the maximum among the
                % minimum values of Vs in the kl-th layer and the minimum
                % values of the shallower layers
                for kl=2:NC
                    var(kl,7)=max(var(kl,7),var(kl-1,7));
                end
                var(:,9)=var(:,8)-var(:,7);
                if find(var(:,9)<0)
                    listerr=find(var(:,9)<0);
                    for kl=1:length(listerr)
                        cap1{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    id=1;
                end
            else
                if find(var(:,9)<0)
                    listerr=find(var(:,9)<0);
                    for kl=1:length(listerr)
                        cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    id=2;
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ZLVp==0
                var(:,4)=max(var(:,4),var(:,7).*sqrt((2*(1-var(:,13)))./(1-(2*var(:,13))))); % We increase minimum Vp if is smaller than the result from minumum poi and minimum Vs
                var(:,5)=min(var(:,5),var(:,8).*sqrt((2*(1-var(:,14)))./(1-(2*var(:,14))))); % We decrease maximum Vp if is greater than the result from maximum poi and maximum Vs
                var(:,6)=var(:,5)-var(:,4);
                % Replacing var(kl,5) (max Vp) with the minimum among the
                % maximum value of Vp in the kl-th layer and the maximun
                % value in the shalower layers
                for kl=NC-1:-1:1
                    var(kl,5)=min(var(kl,5),var(kl+1,5));
                end
                % Replacing var(kl,4) (min Vp) with the maximum among the
                % minimum values of Vp in the kl-th layer and the minimum
                % value in the shallower layers
                for kl=2:NC
                    var(kl,4)=max(var(kl,4),var(kl-1,4));
                end
                var(:,6)=var(:,5)-var(:,4);
                if find(var(:,6)<0)
                    listerr=find(var(:,6)<0);
                    for kl=1:length(listerr)
                        cap1{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    id=1;
                end
            else
                if find(var(:,6)<0)
                    listerr=find(var(:,6)<0);
                    for kl=1:length(listerr)
                        cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    id=2;
                end
            end
        end
        if ZLVp~=0
            if id==1
                errordlg(['Vs values could not be generated for layer ',cap1,' Low velocity zones are forbidden.']);
            elseif id==2
                errordlg(['Range of Vp, Vs and Poisson ratio are incompatible for layer ',cap2]);
            end
        elseif ZLVp==0
            if id==1
                errordlg(['Vp values could not be generated for layer ',cap1,' Low velocity zones are forbidden.']);
            elseif id==2
                errordlg(['Range of Vp, Vs and Poisson ratio are incompatible for layer ',cap2]);
            end
        end
        [dfs]=find(var(:,:)~=varoriginal(:,:));
        if ~isempty(dfs)  && id==0
            if ZLVp~=0
                dia=sprintf('Some values in Vs intervals are incompatible with other constraints!\n\n Upgrade ranges?');
            elseif ZLVp==0
                dia=sprintf('Some values in Vp intervals are incompatible with other constraints!\n\n Upgrade ranges?');
            end
            choice = questdlg(dia, ...
                'Bad condition', ...
                'Update ranges','Cancel','Update ranges');
            % Handle response
            switch choice
                case 'Update ranges'
                    [data(:,5:6),data(:,3:4)]=deal((ceil(var(:,7:8))),(ceil(var(:,4:5))));
                    [handles.TableParameters.Data,fow]=deal(num2cell(data),1);
                    LoadV('var',var,'data',data);
                case 'Cancel'
                    warndlg('Change parameters','!! Warning !!')
                    fow=1;
            end
        end
        [ini1]=GetV('ini1');% Comprobacion del modelo inicial versus parametros de la tabla
        if fow==0 && ~isempty(ini1)
            if find(var(:,1)>ini1(:,1)|var(:,2)<ini1(:,1),1)
                listerr=find(var(:,1)>ini1(:,1)|var(:,2)<ini1(:,1));
                for kl=1:length(listerr)
                    cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                end
                errordlg(['Thickness out of range for layer(s) ',cap2]);
                fow=1;
            elseif find(var(:,4)>ini1(:,2)|var(:,5)<ini1(:,2),1)
                listerr=find(var(:,4)>ini1(:,2)|var(:,5)<ini1(:,2));
                for kl=1:length(listerr)
                    cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                end
                errordlg(['Vp out of range for layer(s) ',cap2]);
                fow=1;
            elseif find(var(:,7)>ini1(:,3)|var(:,8)<ini1(:,3),1)
                listerr=find(var(:,7)>ini1(:,3)|var(:,8)<ini1(:,3));
                for kl=1:length(listerr)
                    cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                end
                errordlg(['Vs out of range for layer(s) ',cap2]);
                fow=1;
            elseif find(var(:,10)>ini1(:,4)|var(:,11)<ini1(:,4),1)
                listerr=find(var(:,10)>ini1(:,4)|var(:,11)<ini1(:,4));
                for kl=1:length(listerr)
                    cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                end
                errordlg(['Density out of range for layer(s) ',cap2]);
                fow=1;
            elseif ZLVp==0 && ~isempty(find(diff(ini1(:,2))<0,1)),
                errordlg(sprintf('Low velocity zones in Vp are forbidden.\n Incompatible initial model.'));
                fow=1;
            elseif ZLVs==0 && ~isempty(find(diff(ini1(:,3))<0,1)),
                errordlg(sprintf('Low velocity zones in Vs are forbidden.\n Incompatible initial model.'));
                fow=1;
            elseif HHS~=0 && ~isempty(find(ini1(end,3)-ini1(:,3)<0,1)),
                errordlg(sprintf('Vs is not maximum for the halfspace in the initial model.\n Check model or settings.'));
                fow=1;
            else
                poi_ini1=(0.5-(ini1(:,3)./ini1(:,2)).^2)./(1-(ini1(:,3)./ini1(:,2)).^2);
                if find(var(:,13)>poi_ini1 | var(:,14)<poi_ini1,1)
                    listerr=find(var(:,13)>poi_ini1 | var(:,14)<poi_ini1);%Poisson
                    for kl=1:length(listerr)
                        cap2{kl}=[' ' num2str(listerr(kl)) ' '];
                    end
                    errordlg(['Poisson`s ratio out of range for layer(s) ',cap2]);
                    fow=1;
                end
            end
        end
        %% fow==0 && id==0 Son consistentes los parametros
        
        if fow==0 && id==0
            % Cleaning memory
            evalin('base','clear(''asens'')');  evalin('base','clear(''ER'')');  evalin('base','clear(''DATAERROR'')');
            evalin('base','clear(''HVALL'')');  evalin('base','clear(''aopt'')');  evalin('base','clear(''MODALL'')');
            evalin('base','clear(''CTALL'')');  evalin('base','clear(''ERROR'')');  evalin('base','clear(''L2ALL'')');
            evalin('base','clear(''L2sens'')');
            % GUI edition
            [handles.edit1.Enable,handles.TDC.Enable,handles.THV.Enable,handles.FHV.Enable,...
                handles.THVDC.Enable,handles.PS.Enable,handles.Models.Enable,handles.Inversion.Enable,...
                handles.LP.Enable,handles.SPara.Enable,handles.smodel.Visible,handles.START.Visible,...
                handles.Report.Enable,handles.Menu.Enable]=deal('off');
            [handles.ploterror.Visible,handles.STOP.Visible]=deal('on');
            %%Inicializacion de variables
            [ER.DCALL, ER.HVALL,ER.asens,ER.L2sens,ER.MODALL,ER.L2ALL,...
                ER.err,ER.err1,ER.temer,ER.kk]=deal([]);
            [handles.STOP.UserData,Cont.iter,Cont.valid,Cont.all]=deal(0);
            
            %% Inversion of HV
            HVDCINV(handles,ER,Cont);
            
            %% GUI edition
            [handles.edit1.Enable,handles.TDC.Enable,handles.THV.Enable,handles.FHV.Enable,...
                handles.THVDC.Enable,handles.PS.Enable,handles.Models.Enable,handles.Inversion.Enable,...
                handles.LP.Enable,handles.SPara.Enable,handles.Menu.Enable,handles.START.Visible ]=deal('on');
            [handles.STOP.Visible,handles.MSA.Checked,handles.SA.Checked,...
                handles.MW.Checked,handles.LS.Checked,handles.LF.Checked]=deal('off');
            Tinv=GetV('Tinv');
            if Tinv==1
                handles.MW.Checked='on';
            elseif Tinv==2
                handles.SA.Checked='on';
            elseif Tinv==3
                handles.MSA.Checked='on';
            elseif Tinv==4
                handles.LS.Checked='on';
            elseif Tinv==5
                handles.LF.Checked='on';
            end
            if ispc
                delete('etc\*.txt');
            elseif ismac
                delete('etc/*.txt');
            elseif isunix
                delete('etc/*.txt');
            end
        end
    else
        errordlg(sprintf('Some values are NAN\n\n Please, change parameters'),'Bad Condition');
    end
else
    errordlg(sprintf('The Poisson ratio values are larger than 0.5\n\n Please, change parameters'),'Bad Condition');
end

%% Stop Inversion
function STOP_Callback(~, ~, handles)
% Funcion for stopping the program during inversion
choice = questdlg('Stop Inversion?','Finished Inversion','Yes','No','Yes');
switch choice
    case 'Yes'
        %Editing GUI
        set(handles.STOP,'Userdata',1);
        set(handles.smodel, 'Visible', 'on')
        set(handles.STOP, 'Visible', 'off');
        set(handles.START, 'Visible', 'on');
        msgbox('Finished program. ','H/V Inversion');
    case 'No'
        handles.STOP.UserData=0;
    otherwise
        handles.STOP.UserData=0;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  "Model Parameters" Structure %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TableParameters_CellEditCallback(hObject, ~, handles)
% Obtaining model parameters ranges from the GUI
data =cell2mat( hObject.Data);
%Fixing halfspace thickness to zero (interpreted as infinity)
if data(end,1)~=0
    data(end,1)=0;
elseif data(end,2)~=0
    data(end,2)=0;
end
%Assigning variables to the Workspace
NC=GetV('NC');
var=zeros(NC,15);
[var(:,1:3:13),var(:,2:3:14)]=deal(data(:,1:2:9),data(:,2:2:10));% Minima to matrix var% Maxima to matrix var
var(:,3:3:15)=var(:,2:3:14)-var(:,1:3:13);
limdd=sum(var(:,2));
LoadV('var',var,'data',data,'limdd',limdd);
handles.TableParameters.Data=num2cell(data);

function LP_Callback(hObject, ~, handles)
%% Import model parameters from *.para
[FileName, Path]=uigetfile({'*.para'},'Open Parameters');
if (FileName)~=0
    fid=fopen([Path FileName]);
    datos=textscan(fid,' %f %f %f %f %f %f %f %f %f %f ','commentstyle','#');
    fclose(fid);
    NC=datos{1}(1);
    if isnan(NC)
        hObject.String=3;
        NC = str2double(hObject.String');
        errordlg('Input must be a number','Error');
    elseif NC<2 || not(mod(NC,1))==0;
        hObject.String= 3;
        NC = str2double(hObject.String);
        errordlg('The number of layers should be greater or equal to 2 and integer','Error');
    else
        [handles.TableParameters.Enable,handles.START.Visible]=deal('on');
        handles.edit1.String=NC;
    end
    data=zeros(NC,10);
    for i=1:10
        for j=1:NC
            data(j,i)= datos{:,i}(j+1);
        end
    end
    var=zeros(NC,15);
    [var(:,1:3:13),var(:,2:3:14)]=deal(data(:,1:2:9),data(:,2:2:10));% Minima to matrix var% Maxima to matrix var
    var(:,3:3:15)=var(:,2:3:14)-var(:,1:3:13);
    [limdd,handles.TableParameters.Data]=deal(sum(var(:,2)),num2cell(data));
    LoadV('limdd',limdd,'var',var,'NC',NC,'data',data);
    [handles.START.Visible,handles.TableParameters.Enable,...
        handles.SPara.Enable,handles.LV.Enable,handles.T0.Enable,...
        handles.Redu.Enable,handles.Inversiont.Enable,handles.edit1.Enable]=deal('on');
end

function SPara_Callback(~, ~, ~)
%% Export model parameters to *.txt
[data1,NC]=GetV('data','NC');
[file,path] = uiputfile('Parameters.para','Save file');
if file~=0
    tex=[path,file];
    dlmwrite(tex,NC);
    dlmwrite(tex, data1 ,'-append', 'delimiter', ' ','precision',6);
end


function NIP_Callback(~, ~, handles)
%% Ask for the nummer of models in the initial population
INPO=GetV('INPO');
if INPO==0
    defaultanswer={num2str(50)};
else
    defaultanswer={num2str(INPO)};
end
NUM=inputdlg({sprintf('Number of models in initial population(>0)\n\n(Introduce 0 to set an initial model)\n')},' ',1,defaultanswer);
if  ~isempty(NUM)
    if ~isnan (str2double(NUM{1}));
        INPO=str2double(NUM{1});
        LoadV('ini1',[],'INPO',INPO);
    end
else
    LoadV('INPO',INPO);
end
if INPO==0
    LoadV('handles',handles);
    initialmodel;
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Funtions for the  Menu %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function PS_Callback(~, ~, ~)
%% Call function which assess sensitivity to parameters
[a,L2]=GetV('MODALL','L2ALL');
sensitiv(a,L2);

function sal_Callback(~, ~, ~)
%% "Exit" (quit program)
choice = questdlg('Quitting?','Are you Sure?','Yes','No','No');
switch choice
    case 'Yes'
        delete('etc\*.txt');
        closereq
    otherwise
        return
end

function Models_Callback(~, ~, ~)
%% New GUI in "Models" menu. Function to show best fitting models
MODELS

%% Two Functions in "Settings" Menu that enable kernels for parallel inversion
function Poff_Callback(~, ~, handles)
open=GetV('open');
if open==1
    [handles.Poff.Checked,handles.Pon.Checked]=deal('on','off');
    mss=msgbox({'Parallelization Settings','Please Wait'},'H/V');
    delete(findobj(mss,'string','OK'));
    delete(findobj(mss,'style','frame'));
    [~]=parallel(1000);
    close(mss);
   GetV('open',0);
    handles.textNumwork.Visible='off';
end

function Pon_Callback(~, ~, handles)
[open,pop]=GetV('open','NUMW');
if open==1
    mss=msgbox({'Parallelization Settings','Please Wait'},'H/V');
    delete(findobj(mss,'string','OK'));
    delete(findobj(mss,'style','frame'));
    delete(gcp('nocreate'))
    close(mss);
end
NUMW=str2double((inputdlg({'Enter number of Workers (Labs): (>=2 & <512)'},'',1,{num2str(pop)})));
if ~isnan(NUMW) || NUMW>=2
    handles.Pon.Checked='on';
    [handles.Poff.Checked,handles.edit1.Enable,handles.TDC.Enable,handles.THV.Enable,...
     handles.FHV.Enable,handles.THVDC.Enable,handles.PS.Enable,handles.Models.Enable,...
     handles.Inversion.Enable,handles.LP.Enable,handles.SPara.Enable,handles.START.Enable,...
     handles.STOP.Enable,handles.Menu.Enable]=deal('off');

    mss=msgbox({sprintf('The Parallelization Settings will take a few seconds.\n\nPlease Wait')},'HV-inv');
    delete(findobj(mss,'string','OK'));
    delete(findobj(mss,'style','frame'));
    set(mss, 'CloseRequestFcn','')
    parallel(NUMW);
    delete(mss);
    [handles.textNumwork.String,handles.edit1.Enable]=...
        deal(['Connection with',' ',num2str(NUMW),' ', 'Labs.'],'off');
 
    set(handles.textNumwork,'Visible','on');
    set(handles.TDC, 'Enable', 'on');
    set(handles.THV, 'Enable', 'on');
    set(handles.FHV, 'Enable', 'on');
    set(handles.THVDC, 'Enable', 'on');
    set(handles.PS, 'Enable', 'on');
    set(handles.Models, 'Enable', 'on');
    set(handles.Inversion, 'Enable', 'on');
    set(handles.LP, 'Enable', 'on');
    set(handles.SPara, 'Enable', 'on');
    set(handles.START, 'Enable', 'on');
    set(handles.STOP, 'Enable', 'on');
    set(handles.Menu, 'Enable', 'on');
    
    LoadV('NUMW',NUMW,'open',1);
end

function T0free_Callback(~, ~,  handles)
%% Function for user defined "Initial Temperature" menu
AP=str2double((inputdlg({'Temperature'},'T0',1,{'10'})));
LoadV('AP',AP,'T0',2)
[handles.T0free.Checked,handles.T0misfit.Checked]=deal('on','off');

function T0misfit_Callback(~, ~,handles)
%% Function for "Initial Temperature" menu (alternative input method)
T0=1;
defaultanswer={'0.1','0.5'};
prompt={'Relative misfit increment (>0)',' Probability of acceptance (<1)'};
name='T0';
numlines=1;
NUM=inputdlg(prompt,name,numlines,defaultanswer);
if ~isempty(NUM)
    EA=str2double(NUM{1});
    AP=str2double(NUM{2});
    assignin('base','EA',EA)
    assignin('base','AP',AP)
    set(handles.T0free,'Checked','off')
    set(handles.T0misfit,'Checked','on')
    assignin('base','T0',T0)
end

function Redu_Callback(~, ~, ~)
%% Function for "Cooling Schedule" menu
defaultanswer={'0.9'};
prompt={'Temperature ratio (Ti+1/Ti)'};
name='Cooling Schedule';
numlines=1;
RT=str2double(inputdlg(prompt,name,numlines,defaultanswer));
assignin('base','RT',RT)

function MW_Callback(~, ~, handles)
%% Function for "Algorithm details" menu
defaultanswer={'100','10'};
prompt={'Number of iterations','Perturbation range (%)'};
name=' ';
numlines=1;
NUM=((inputdlg(prompt,name,numlines,defaultanswer)));
if ~isempty(NUM)
    A=evalin('base','A');
    set(handles.MSA,'Checked','off')
    set(handles.SA,'Checked','off')
    set(handles.MW,'Checked','on')
    set(handles.LS,'Checked','off')
    set(handles.LF,'Checked','off')
    set(handles.Redu,'Enable','off')
    set(handles.T0,'Enable','off')
    Tinv=1;
    A.Nsubiter=0;
    A.Niter=0;
    A.Niterf=str2double((NUM{1}));
    A.NN=str2double((NUM{2}));
    assignin('base','A',A)
    assignin('base','Tinv',Tinv)
end

function SA_Callback(~, ~, handles)
%% Function for Simulated Annealing (SA) menu
defaultanswer={'100','0','100','10'};
prompt={['M',char(225),'rkov`s chain length'],'Number of temperatures',['Last M',char(225),'rkov`s chain length'],'Perturbation range (%)'};
name=' ';
numlines=1;
NUM=inputdlg(prompt,name,numlines,defaultanswer);
if ~isempty(NUM)
    A=evalin('base','A');
    set(handles.MSA,'Checked','off')
    set(handles.SA,'Checked','on')
    set(handles.MW,'Checked','off')
    set(handles.Redu,'Enable','on')
    set(handles.T0,'Enable','on')
    set(handles.LS,'Checked','off')
    set(handles.LF,'Checked','off')
    Tinv=2;
    A.Nsubiter=str2double((NUM{2}));
    A.Niter=str2double((NUM{1}));
    A.Niterf=str2double((NUM{3}));
    A.NN=str2double((NUM{4}));
    assignin('base','A',A)
    assignin('base','Tinv',Tinv)
end


function MSA_Callback(~, ~, handles)
%% Function for Modified Simulated Annealing (MSA) menu
defaultanswer={'100','0','100','10'};
prompt={'Number of iterations','Number of Reheatings','Number of Last iterations','Perturbation range (%)'};
name=' ';
numlines=1;
NUM=((inputdlg(prompt,name,numlines,defaultanswer)));
if ~isempty(NUM)
    A=evalin('base','A');
    set(handles.MSA,'Checked','on')
    set(handles.SA,'Checked','off')
    set(handles.MW,'Checked','off')
    set(handles.Redu,'Enable','on')
    set(handles.T0,'Enable','on')
    set(handles.LS,'Checked','off')
    set(handles.LF,'Checked','off')
    Tinv=3;
    A.Nsubiter=str2double((NUM{2}));
    A.Niter=str2double((NUM{1}));
    A.Niterf=str2double((NUM{3}));
    A.NN=str2double((NUM{4}));
    assignin('base','A',A);
    assignin('base','Tinv',Tinv)
end
open=evalin('base','open');
if open==0
    choice = questdlg(sprintf('The inversion method is designed to run in parallel.\n\nTo use this method in parallel?'));
    if strcmp('Yes',choice);
        
        pop=evalin('base','NUMW');
        defaultanswer={num2str(pop)};
        prompt={sprintf('Enter number of Workers (Labs): (>=2 & <512)')};
        name=' ';
        numlines=1;
        NUMW=str2double(cell2mat(inputdlg(prompt,name,numlines,defaultanswer)));
        if isnan(NUMW) || NUMW<2
        else
            set(handles.Poff,'Checked','off');
            set(handles.Pon,'Checked','on');
            mss=msgbox({sprintf('The Parallelization Settings will take a few seconds.\n\nPlease Wait')},'HV-inv');
            delete(findobj(mss,'string','OK'));
            delete(findobj(mss,'style','frame'));
            set(mss, 'CloseRequestFcn','')
            set(handles.edit1, 'Enable', 'off');
            set(handles.TDC, 'Enable', 'off');
            set(handles.THV, 'Enable', 'off');
            set(handles.FHV, 'Enable', 'off');
            set(handles.THVDC, 'Enable', 'off');
            set(handles.PS, 'Enable', 'off');
            set(handles.Models, 'Enable', 'off');
            set(handles.Inversion, 'Enable', 'off');
            set(handles.LP, 'Enable', 'off');
            set(handles.SPara, 'Enable', 'off');
            set(handles.START, 'Enable', 'off');
            set(handles.STOP, 'Enable', 'off');
            set(handles.Menu, 'Enable', 'off');
            parallel(NUMW);
            delete(mss);
            open=1;
            assignin('base','NUMW',NUMW);
            assignin('base','open',open);
            set(handles.textNumwork,'Visible','on');
            set(handles.textNumwork,'String',['Connection with',' ',num2str(NUMW),' ', 'Labs.']);
            set(handles.edit1, 'Enable', 'off');
            set(handles.TDC, 'Enable', 'on');
            set(handles.THV, 'Enable', 'on');
            set(handles.FHV, 'Enable', 'on');
            set(handles.THVDC, 'Enable', 'on');
            set(handles.PS, 'Enable', 'on');
            set(handles.Models, 'Enable', 'on');
            set(handles.Inversion, 'Enable', 'on');
            set(handles.LP, 'Enable', 'on');
            set(handles.SPara, 'Enable', 'on');
            set(handles.START, 'Enable', 'on');
            set(handles.STOP, 'Enable', 'on');
            set(handles.Menu, 'Enable', 'on');
        end
    end
end

function LS_Callback(~, ~, handles)
%% Function for Simplex Downhill menu
NUM=str2double(inputdlg({'Maximum number of iterations allowed','Maximum number of evaluations allowed',...
    'Termination tolerance on the function value (TolFun)','Termination tolerance on parameters (TolX)'},...
    ' ',1,{'100','200','1e-4','1e-4'}));
if ~isempty(NUM)
    [handles.MSA.Checked,handles.SA.Checked,handles.MW.Checked,...
        handles.LF.Checked,handles.Redu.Enable,handles.T0.Enable]=deal('off');
    [A.Niter,A.Nsubiter,A.Niterf,A.NN,handles.LS.Checked]=deal(NUM(1),NUM(2),NUM(3),NUM(4),'on');
    LoadV('Tinv',4,'A',A)
end

%% Function for Interior point (IP) menu
function LF_Callback(~, ~, handles)
NUM=str2double(inputdlg({'Maximum number of iterations allowed','Maximum number of evaluations allowed',...
    'Termination tolerance on the function value (TolFun)','Termination tolerance on parameters (TolX)'},...
    ' ',1,{'100','200','1e-4','1e-4'}));
if ~isempty(NUM)
    [handles.MSA.Checked,handles.SA.Checked,handles.MW.Checked,...
       handles.LS.Checked,handles.Redu.Enable,handles.T0.Enable]=deal('off');
    [A.Niter,A.Nsubiter,A.Niterf,A.NN,handles.LF.Checked]=deal(NUM(1),NUM(2),NUM(3),NUM(4),'on');
    LoadV('Tinv',5,'A',A)
end

%% Function for the "Wave Parameters" menu
function WPara_Callback(~, ~, ~)
NUM=str2double(inputdlg({'Rayleigh waves modes','Love wave modes','Minimum number of integration points for BW' ,'Maximum Number of integration points for BW','Regularization factor'},...
    'Waves Parameters',1,{'5','5','1000','2000','0.001'}));
if ~isempty(NUM)
    [NR,NL,dk,mdk,apsv]=deal(NUM(1),NUM(2),NUM(3),NUM(4),NUM(5));
    if NR<0 || NL <0
        errordlg('The minimum number of waves modes should be greater to zero','Error');
    end
    if dk>mdk
        errordlg('The maximun number of integration should be greater or equal to minimum number of integration ','Error');
        dk=mdk;
    end
    LoadV('ramasR',NR,'ramasL',NL,'kb_by_dk',dk,'mkb_by_dk',mdk,'apsv',apsv)
end

%% Function for quitting the program
function figure1_CloseRequestFcn(~, ~, ~)
choice = questdlg('Quitting?','Are you Sure?','Yes','No','No');
switch choice
    case 'Yes'
        delete('etc\*.txt');
        closereq
    otherwise
        return
end

%% Three functions for the "L_V Zones" Menu
%% Function for Low S-wave velocity zones
function LZVS_ButtonDownFcn(~, ~, handles)
ZLVp=evalin('base','ZLVp');
if strcmp(get(handles.LZVS,'Checked'),'off')
     [handles.LZVS.Checked,ZLVs]=deal('on',1);
else
        [handles.LZVS.Checked,ZLVs]=deal('off',0);
    %     HHS=1; set(handles.HHS,'Checked','on')
end
if ZLVp==0 && ZLVs==0
    [handles.LZVP.Checked,ZLVp]=deal('on',1);
    %     HHS=0; set(handles.HHS,'Checked','off')
end
LoadV('ZLVs',ZLVs,'ZLVp',ZLVp);
% assignin('base','HHS',HHS);

%% Function for Low P-wave velocity zones
function LZVP_ButtonDownFcn(~, ~, handles)
ZLVs=GetV('ZLVs');
if strcmp(get(handles.LZVP,'Checked'),'off')
    [handles.LZVP.Checked,ZLVp]=deal('on',1);
else
    [handles.LZVP.Checked,ZLVp]=deal('off',0);
end
if ZLVp==0 && ZLVs==0
    [handles.LZVS.Checked,ZLVs]=deal('on',1);
end
LoadV('ZLVs',ZLVs,'ZLVp',ZLVp)

%% Function for Hard half space
function HHS_ButtonDownFcn(~, ~, handles)
if strcmp(get(handles.HHS,'Checked'),'off')
    [handles.HHS.Checked,handles.LZVS.Checked,HHS,ZLVs]=...
        deal('on','off',1,0);
else
     [handles.HHS.Checked,handles.LZVS.Checked,HHS,ZLVs]=...
        deal('off','on',0,1);
end
LoadV('HHS',HHS,'ZLVs',ZLVs);

function varargout = HVTI_OutputFcn(~, ~, handles)
%% Function for "Help-About" menu
varargout{1} = handles.output;
msgbox(sprintf(['HVInv 1.31 Beta Project 2014-2016\n\n',...
    'This program is free software: you can redistribute it and/or modify\n',...
    'it under the terms of the GNU General Public License version 3 as\n',...
    'published by the Free Software Foundation. For details see the copy\n',...
    'included in this pakage.\n\n',...
    'This program is distributed in the hope that it will be useful,\n',...
    'but WITHOUT ANY WARRANTY; without even the implied warranty of\n',...
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\nSee the ',...
    'GNU General Public License for more details.\n\n',...
    'You should have received a copy of the GNU General Public License\n',...
    'along with this program.  If not, see <http://www.gnu.org/licenses/>..\n\n',...
    'Contributors:\n\n',...
    'Antonio Garc',char(237),'a-Jerez (agarcia-jerez@ual.es)\n',...
    'Jos',char(233), ' Pi',char(241),'a-Flores (ead2009@hotmail.com)\n',...
    'Mathieu Perton (mathieu.perton@gmail.com)\n',...
    'Francisco J. S',char(225),'nchez-Sesma (sesma@unam.mx)\n',...
    'Francisco Luz',char(243),'n (fluzon@ual.es)\n',...
    'Marc Wathelet\n\n',...
    'Forward calculation based on original ideas of:\n\n',...
    'F. S',char(225),'nchez-Sesma, A. Garc',char(237),'a-Jerez and M. Perton\n\n',...
    'Support for GUI and inverse problem:\n\nJos',char(233), ' Pi',char(241),'a-Flores\n\n',...
    'Support for forward calculations:\n\nA. Garc',char(237),'a-Jerez\n\n',...
    'More information at http://www.ual.es/GruposInv/hv-inv/\n\n'...
    'MATLAB' char(174),',', 'Version: 8.5.0.197613 (R2015a) &  9.0.0 341360 (R2016a)\n'...
    'License Numbers: 125381, 40274058, 40462394 \n\n\n'...
    char(169),' Jos',char(233), ' Pi',char(241),'a-Flores & Antonio Garc',char(237),'a-Jerez, 2014-2016\n']),'About');


function About_Callback(~, ~, ~)
msgbox(sprintf(['HVInv 1.31 Beta Project 2014-2016\n\n',...
    'This program is free software: you can redistribute it and/or modify\n',...
    'it under the terms of the GNU General Public License version 3 as\n',...
    'published by the Free Software Foundation. For details see the copy\n',...
    'included in this pakage.\n\n',...
    'This program is distributed in the hope that it will be useful,\n',...
    'but WITHOUT ANY WARRANTY; without even the implied warranty of\n',...
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\nSee the ',...
    'GNU General Public License for more details.\n\n',...
    'You should have received a copy of the GNU General Public License\n',...
    'along with this program.  If not, see <http://www.gnu.org/licenses/>..\n\n',...
    'Contributors:\n\n',...
    'Antonio Garc',char(237),'a-Jerez (agarcia-jerez@ual.es)\n',...
    'Jos',char(233), ' Pi',char(241),'a-Flores (ead2009@hotmail.com)\n',...
    'Mathieu Perton (mathieu.perton@gmail.com)\n',...
    'Francisco J. S',char(225),'nchez-Sesma (sesma@unam.mx)\n',...
    'Francisco Luz',char(243),'n (fluzon@ual.es)\n',...
    'Marc Wathelet\n\n',...
    'Forward calculation based on original ideas of:\n\n',...
    'F. S',char(225),'nchez-Sesma, A. Garc',char(237),'a-Jerez and M. Perton\n\n',...
    'Support for GUI and inverse problem:\n\nJos',char(233), ' Pi',char(241),'a-Flores\n\n',...
    'Support for forward calculations:\n\nA. Garc',char(237),'a-Jerez\n\n',...
    'More information at http://www.ual.es/GruposInv/hv-inv/\n\n'...
    'MATLAB' char(174),',', 'Version: 8.5.0.197613 (R2015a) &  9.0.0 341360 (R2016a)\n'...
    'License Numbers: 125381, 40274058, 40462394 \n\n\n'...
    char(169),' Jos',char(233), ' Pi',char(241),'a-Flores & Antonio Garc',char(237),'a-Jerez, 2014-2016\n']),'About');

%% Open User manual
function Manual_Callback(~, ~, ~)
open('./Document/UserManual.html')

function smodel_Callback(~, ~, ~)
%% Save optimum model in *.txt format
aopt=evalin('base', 'aopt');
[file,path] = uiputfile('Model.txt','Save file name');
if ~isempty(file)
    tex=[path,file];
    fileID = fopen(tex,'w');
    for row = 1:length(aopt(:,1))+1
        if row==1
            fprintf(fileID,'%s\r' ,num2str(length(aopt(:,1))));
        else
            fprintf(fileID,'%s\r' ,num2str(aopt(row-1,:)));
        end
    end
    fclose(fileID);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Additional GUI Functions %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NIP_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



