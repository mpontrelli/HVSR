%     Copyright (C) 2014,2016 José Piña-Flores, Antonio García-Jerez.
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License version 3 as
%     published by the Free Software Foundation.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more deta
function varargout = FHV(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @FHV_OpeningFcn, ...
    'gui_OutputFcn',  @FHV_OutputFcn, ...
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
function FHV_OpeningFcn(hObject, ~, handles, varargin)
global ground_par keep_poi LockAxis change1thickness keep_hv %parametros globales


% GUI default Variables
handles.output = hObject;% Control de los elementos de la GUI
guidata(hObject, handles);% Control de los elementos de la GUI
[handles.uitable1.Enable,handles.fmin.Enable,handles.fmax.Enable,handles.L.Enable,...
    handles.R.Enable,handles.bwi.Enable,handles.save.Enable,handles.samples.Enable,...
    handles.pushbutton1.Enable,handles.LMHV.Enable]=deal('off');% Control de los elementos de la GUI
set(handles.uipanel12.Children(:),'Enable','off');% Control de los elementos de la GUI

% Default settings for calling the forward H/V program:
% Frequency,samples,Modes SW and BW
[handles.fmin.String,handles.fmax.String,handles.samples.String,...
    handles.bwi.String,handles.L.String,handles.R.String]=deal(0.2,10,100,500,10,10);
% Global variables (Elastic parameter)
[ground_par,change1thickness,LockAxis,keep_poi,keep_hv]=deal(3,0,0,false,false);
% Identificator of graphical panels in GUI
% Panel for HV
spHV=subplot(1,1,1,'Parent',handles.uipanel11);% subplot for HVs
spHV.NextPlot='add';
ThHVid=plot(spHV,nan,nan,'b','LineWidth',3);% Theoretical HV ID
ExHVid=errorbar(spHV,nan,nan,nan,'--r','LineWidth',2);% Experimental HV
ThHVBM=plot(spHV,nan,nan,'-k','LineWidth',3);% Theoretical HV of loaded (back) model
% Activate Zoom and Pan tool into HV and Profile grafic
ThHVid.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));
ExHVid.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));
ThHVBM.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));
ZoomAction=zoom;
ZoomAction.UIContextMenu=uicontextmenu;
ZoomAction.UIContextMenu.Callback=@(hObject,eventdata)FHV('DesactiveTool',hObject,eventdata,guidata(hObject));
set(get(ExHVid,'Children'),{'LineWidth'},{2; 0.5})
% Graphical Profile
spMDL=subplot(1,1,1,'Parent',handles.uipanel14);% subplot for the ground model
spMDL.NextPlot='add';
auxMDLid=plot(spMDL,nan,nan,'.-b','LineWidth',3);% Auxiliary profile
% Asignacion de la funcion para calcular el HV con un click del mouse
auxMDLid.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));
% Menu context for graphic
cF = uicontextmenu; %contexmenu for HV
cP=uicontextmenu;% contexmenu for model
[ExHVLegend,ThHVBMLegend,ClearHVDataOption]=deal(' '); % legends
[spHV.XLabel.String,spHV.YLabel.String,spMDL.XLabel.String,spMDL.YLabel.String]=deal...
    ('Frequency [Hz]','Amplitude (${H} \over {V}$)','Velocity Vs [m/s]','Depth [m]');% Labels
[spHV.XLabel.Interpreter,spHV.YLabel.Interpreter,spMDL.XLabel.Interpreter,spMDL.YLabel.Interpreter]=deal('latex');
[spHV.FontUnits,spMDL.FontUnits]=deal('normalized');
% Load default variables
LoadV('fminFHV',0.2,'fmaxFHV',10,'nFHV',100,'nksBW',500,'LFHV',10,'RFHV',10,'IsLogF',1,'apsvFHV',0.001,'ShowAuxMDL',0,...
    'bmodel',nan,'limdd',nan,'HVFHV',nan,'FFHV',nan,'HVFHVD',nan,'FFHVD',nan,'spHV',spHV,'spMDL',spMDL,'ThHVid',ThHVid,'ThHVBM',ThHVBM,...
    'auxMDLid',auxMDLid,'ExHVid',ExHVid,'ExHVLegend',ExHVLegend,'ThHVBMLegend',ThHVBMLegend','cF',cF,...
    'ClearHVDataOption',ClearHVDataOption,'ZoomAction',ZoomAction,'ZoomOn',0);
[spHV.UIContextMenu,spMDL.UIContextMenu]=deal(cF,cP);
% Legend for Menu context for HV
uimenu(cF,'Label','XScale: Log','Callback',@setlinestyle);
uimenu(cF,'Label','XScale: Linear','Callback',@setlinestyle);
uimenu(cF,'Label','YScale: Log','Callback',@setlinestyle);
uimenu(cF,'Label','YScale: Linear','Callback',@setlinestyle);
% % Legend for Menu context for model
uimenu(cP,'Label','Lock X Axis','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Lock Y Axis','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Keep Poisson ratio','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Change one thickness','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Move one interface','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Vp Profile','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Vs Profile','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Density Profile','Callback',@HandleMouseMenu);
uimenu(cP,'Label','Keep HV ratio','Callback',@HandleMouseMenu);
LoadV('cP',cP);
% Apariencia del menu context
[cP.Children(3).Checked,cP.Children(5).Checked,cP.Children(4).Separator,...
    cP.Children(6).Separator,cP.Children(1).Separator,cF.Children(2).Separator]=deal('on');
[spMDL.Visible,spHV.Visible]=deal('off');
% Graphic control of the menu GUI
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


%% GUI Function edit minimum frequency by user
function fmin_Callback(hObject, ~, ~)
fminFHV = str2double(hObject.String); % Obtiene el valor de la frecuencia minima introducido por el usuario
if isnan(fminFHV)
    hObject.String=0.2;
    errordlg('Input must be a number','Error');
elseif fminFHV<=0
    hObject.String=0.2;
    errordlg('The number of minimun frequency must be greater than 0','Error');
end
LoadV('fminFHV',fminFHV);% Asignacion del valor al workspace

%% GUI Function edit maximum frequency by user
function fmax_Callback(hObject, ~, ~)
fminFHV=GetV('fminFHV');%Get fmin variable from workspace
fmaxFHV = str2double(hObject.String);%Get fmax variable from gui
if isnan(fmaxFHV)
    [hObject.String,fmaxFHV]=deal(10);
    errordlg('Input must be a number','Error');
elseif fmaxFHV<=fminFHV
    [hObject.String,fmaxFHV]=deal(10);
    errordlg('The number of maximun frequency must be greater than minimum frequency','Error');
end
LoadV('fmaxFHV',fmaxFHV);% Asignacion del valor de frecuencia maxima al workspace

%% GUI Function edit number of samples
function samples_Callback(hObject, ~, ~)
nF = str2double(hObject.String);% Número de muestras
%Condiciones
if isnan(nF)
    [hObject.String,nF]=deal(100);
    errordlg('Input must be a number','Error');
elseif nF<2 || not(mod(nF,1))==0;
    [hObject.String,nF]=deal(100);
    errordlg('The number of samples must be greater than 2 and integer','Error');
end
LoadV('nFHV',nF);% Load variable

%% GUI sampling options
function uipanel12_SelectionChangeFcn(hObject,~, handles)
if hObject==handles.Log%get variable
    IsLogF=1;
else
    IsLogF=0;
end
LoadV('IsLogF',IsLogF);%Load Variable

%% GUI Function edit Body wave integrals by user
function bwi_Callback(hObject, ~, ~)
nksBW = str2double(hObject.String);% Get variable
%Condiciones
if isnan(nksBW)
    [hObject.String,nksBW]=deal(50);
    errordlg('Input must be a number','Error');
elseif nksBW<0 || not(mod(nksBW,1))==0;
    [hObject.String,nksBW]=deal(50);
    errordlg('The number of integration must be positive and integer','Error');
end
LoadV('nksBW',nksBW);% Load Variable

%% GUI Function edit number of Love modes by user
function L_Callback(hObject, ~, ~)
LF = str2double(hObject.String);% Get Variable
if isnan(LF)
    [hObject.String,LF]=deal(10);
    errordlg('Input must be a number','Error');
elseif LF<0  || not(mod(LF,1))==0;
    [hObject.String,LF]=deal(10);
    errordlg('The number of modes must be greater than 0 and integer','Error');
end
LoadV('LFHV',LF);% Load Variable

%% GUI Function edit number of Rayleigh modes by user
function R_Callback(hObject, ~, ~)
RF = str2double(hObject.String);%Get Variable
%Condiciones
if isnan(RF)
    [hObject.String,RF]=deal(10);
    errordlg('Input must be a number','Error');
elseif RF<0  || not(mod(RF,1))==0;
    [hObject.String,RF]=deal(10);
    errordlg('The number of modes must be greater than 0 and integer','Error');
end
LoadV('RFHV',RF);%Load Variable

%% GUI Function save file of HV data
function save_Callback(~, ~, ~)
[HV,FFHV]=GetV( 'HVFHV', 'FFHV');%Load data from theorical HV
aop=[FFHV HV];
[file,path] = uiputfile('HV.txt','Save file name');%Path
if file~=0
    tex=[path,file];
    dlmwrite(tex, aop, 'delimiter', ' ','precision',6);% Save data
end

%% GUI Stabilization factor (slight damping for body wave integrals)
function RFact_Callback(~, ~, ~)
NUM=(inputdlg('Regularization factor','Waves Parameters',1,{'0.01'}));
if ~isempty(NUM)
    apsvFHV=str2double(NUM{1});
else
    apsvFHV=0.01;
end
LoadV('apsvFHV',apsvFHV);% Load Variable

%% GUI Function elastic parameters to default
function edit1_Callback(hObject, ~, handles)
global dataFHV
NCF = str2double(hObject.String);
if isnan(NCF) || isinf(NCF);
    [hObject.String,NCF]=deal(2);
    errordlg('Input must be a number','Error');
elseif NCF<1 || not(mod(NCF,1))==0;
    hObject.String= 2;
    NCF = 2;
    errordlg('The number of layers should be greater or equal to 1 and integer','Error');
end
% Control de GUI
handles.uipanel14.Children.Visible='on';
[handles.uitable1.Enable,handles.pushbutton1.Enable,handles.fmin.Enable,handles.fmax.Enable,...
    handles.samples.Enable,handles.bwi.Enable,handles.L.Enable,handles.R.Enable,handles.SM.Enable,...
    handles.uipanel12.Children(:).Enable]=deal('On');
handles.LMHV.Value=0;
handles.LMHV.Enable= 'Off';
% Valores de la tabla por default
numeleF=zeros(NCF,5);
[numeleF(:,1),numeleF(end,1),numeleF(:,2),numeleF(end,2),...
    numeleF(:,3),numeleF(end,3),numeleF(:,4),numeleF(end,4)] =deal(10,0,500,1000,200,500,1500,2500);%thickness
numeleF(:,5)=(0.5-(numeleF(:,3)./numeleF(:,2)).^2)./(1-(numeleF(:,3)./numeleF(:,2)).^2);%Poisson
%Asignando valores
[dataFHV,handles.uitable1.Data,handles.uitable1.ColumnEditable]=deal(numeleF,num2cell(numeleF),true(1,5));
% Limpiando el modelo backmodel
[ThHVBM,auxMDLid]=GetV('ThHVBM','auxMDLid');
[ThHVBM.XData,ThHVBM.YData]=deal(nan);

if length(auxMDLid.XData)~=1% 1= No exite propiedades de impoly (g) ~=1 exite g
    g=GetV('g'); % Load local copy of g from workspace
    delete((g.Parent.Children(1)));% delete impoly without line
end
% Llamando funcion interactiva para modelo
VisualFit_v8(auxMDLid,handles);
g=GetV('g');
for ii=1:2*NCF
    delete(g.Children(ii).UIContextMenu.Children)
end;% borra el menu del botón derecho
delete(g.Children(end).UIContextMenu.Children)
[g.Parent.Children(end).XData,g.Parent.Children(end).YData]=deal(nan(size(g.Parent.Children(end).YData)));
% [g.Parent.Children(end).Marker,g.Parent.Children(end).LineStyle]=deal('none');
LIM=sum(dataFHV(:,1))*1.25;
handles.uipanel14.Children.YLim=[0 sum(dataFHV(:,1))*1.25 ];
LoadV('bmodel',nan,'limdd',nan,'NCFHV',NCF,'ZoomOn',0,'LIM',LIM)% Load Variables

%% GUI Function edit elastic parameters by user (Table of parameters)
function uitable1_CellEditCallback(hObject, ~, handles)
global ground_par keep_poi LockAxis dataFHV pol_backup
[dataF, dataFa] = deal(cell2mat(hObject.Data));% Obteniendo valores de la tabla
dataFb=dataFHV;% Asignando a variable local
[auxMDLid,ShowAuxMDL,g,bmodel,h]=GetV('auxMDLid','ShowAuxMDL','g','bmodel','h');%perfil
% Formato para pintar modelo
[B,A,R,E]=MODELO(bmodel);
LockAxisBack=LockAxis;%LockAxis=evalin('base','LockAxis');
if isempty(find(dataFa(1:end-1,1)~=dataFb(1:end-1,1),1)),% No thickness changed
    if ~isempty(find(dataF(:,5)~=dataFb(:,5),1)),% Poisson's ratio changed, update Vp
        ground_par=2;% Model caso Vp
        %Apariencia de GUI
        [auxMDLid.Parent.XLabel.String,auxMDLid.Parent.Parent.Title]=deal('Velocity Vp [m/s]','Vp Profile');
        auxMDLid.Parent.UIContextMenu.Children(4).Checked='on';
        [auxMDLid.Parent.UIContextMenu.Children(3).Checked,...
            auxMDLid.Parent.UIContextMenu.Children(2).Checked]=deal('off');
        % Calculo del coeficiente de poisson
        if ismac
            dataF(:,2)=round(real(dataF(:,3).*(sqrt((2*(1-dataF(:,5)))./(1-(2.*dataF(:,5)))))),2);
        else
            dataF(:,2)=roundn(real(dataF(:,3).*(sqrt((2*(1-dataF(:,5)))./(1-(2.*dataF(:,5)))))),-2);
        end
    elseif keep_poi, % A velocity changed and we want to preserve Poisson's ratio
        if ismac
            if ~isempty(find(dataF(:,3)~=dataFb(:,3),1)),% Vs changed, update Vp
                dataF(:,2)=round(real(dataF(:,3).*(sqrt((2*(1-dataF(:,5)))./(1-(2.*dataF(:,5)))))),2);
            elseif ~isempty(find(dataF(:,2)~=dataFb(:,2),1)) % Vp changed, update Vs
                dataF(:,3)=round(real(dataF(:,2).*(sqrt((1-(2.*dataF(:,5)))./(2*(1-dataF(:,5)))))),2);
            end
        else
            if ~isempty(find(dataF(:,3)~=dataFb(:,3),1)),% Vs changed, update Vp
                dataF(:,2)=roundn(real(dataF(:,3).*(sqrt((2*(1-dataF(:,5)))./(1-(2.*dataF(:,5)))))),-2);
            elseif ~isempty(find(dataF(:,2)~=dataFb(:,2),1)) % Vp changed, update Vs
                dataF(:,3)=roundn(real(dataF(:,2).*(sqrt((1-(2.*dataF(:,5)))./(2*(1-dataF(:,5)))))),-2);
            end
        end
    else % A velocity changed and we want to preserve the other velocity
        if ismac
            dataF(:,5)=round((0.5-(dataF(:,3)./dataF(:,2)).^2)./(1-(dataF(:,3)./dataF(:,2)).^2),3);
        else
            dataF(:,5)=roundn((0.5-(dataF(:,3)./dataF(:,2)).^2)./(1-(dataF(:,3)./dataF(:,2)).^2),-3);
        end
    end
end
[dataF(end,1),dataFHV,handles.uitable1.Data]=deal(0,dataF,num2cell(dataF));
% Change the graphic model to show the variable changed in the table
if ~isempty(find(dataFb(:,2)~=dataFa(:,2),1))
    %Since the user changed Vp, select Vp in the graphic model
    %Asignando valores a variables
    [ground_par,auxMDLid.Parent.XLabel.String,auxMDLid.Parent.Parent.Title,auxMDLid.Parent.UIContextMenu.Children(4).Checked,...
        auxMDLid.Parent.UIContextMenu.Children(3).Checked,auxMDLid.Parent.UIContextMenu.Children(2).Checked]=...
        deal(2,'Velocity Vp [m/s]','Vp Profile','on','off','off');
    if ShowAuxMDL
        [g.Parent.Children(end).XData,g.Parent.Children(end).YData]=deal(A,E);
    end
elseif sum(dataFb(:,3)~=dataFa(:,3))
    %The user changed Vs, select Vs in the graphic model
    [ground_par,auxMDLid.Parent.XLabel.String,auxMDLid.Parent.Parent.Title,auxMDLid.Parent.UIContextMenu.Children(4).Checked,...
        auxMDLid.Parent.UIContextMenu.Children(3).Checked,auxMDLid.Parent.UIContextMenu.Children(2).Checked]=...
        deal(3,'Velocity Vs [m/s]','Vs Profile','off','on','off');
    if ShowAuxMDL
        [g.Parent.Children(end).XData,g.Parent.Children(end).YData]=deal(B,E);
    end
elseif sum(dataFb(:,4)~=dataFa(:,4))
    %The user changed density, select density in the graphic model
    [ground_par,auxMDLid.Parent.XLabel.String,auxMDLid.Parent.Parent.Title,auxMDLid.Parent.UIContextMenu.Children(4).Checked,...
        auxMDLid.Parent.UIContextMenu.Children(3).Checked,auxMDLid.Parent.UIContextMenu.Children(2).Checked]=...
        deal(4,'Density [kg/m3]','Density Profile','off','off','on');
    if ShowAuxMDL
        [g.Parent.Children(end).XData,g.Parent.Children(end).YData]=deal(R,E);
    end
end
% % Formato para pintar modelo
[B,A,R,E]=MODELO(dataF(:,1:4));
%Actualizando cambios en la grafica del perfil
E(end)=sum(dataF(:,1))*1.25;
ProfilePolys=[E;A;B;R].';
pol_backup=[ProfilePolys(:,ground_par) ProfilePolys(:,1)];
setConstrainedPosition(h,pol_backup);
LockAxis=LockAxisBack;% Restore LockAxis status

%% GUI Function for start computing HV curve Button
function pushbutton1_Callback(~, ~, handles)
CalHV(handles,1);

%% Function for calculate HV and activate Zoom Pan tools by mouse click into Profile
function ClickAction(~,event,handles)
persistent click
if event.Button==1 %Action with the primary button
    if isempty(click)
        click = 1;
        pause(0.5); %Add a delay to distinguish single click from a double click
        if click == 1
            click = [];
        end
    else
        CalHV(handles,1)% Call Calculate HV Theoric
        click = [];
    end
elseif event.Button==2 % Action with the scroll button
    if isempty(click)
        click = 1;
        pause(0.5); %Add a delay to distinguish single click from a double click
        if click == 1
            % Assignin Uicontextmenu Pan Tool
            ZoomAction=GetV('ZoomAction');
            ZoomAction.Direction='in';
            ZoomAction.Enable='on';
            LoadV('ZoomOn',1);
            click = [];
        end
    else
        if  strcmp(event.Source.Parent.Title,'HV ratio')
            axis(handles.uipanel11.Children(end),'tight');% HV curve axis
            axis(handles.uipanel11.Children(end),'manual');
        else %% Model axis
            axis(handles.uipanel14.Children,'tight');
            LoadV('ZoomOn',0);
            axis(handles.uipanel14.Children,'manual');
        end
        click = [];
    end
end

%% Desactivation Pan and Zoom tool
function DesactiveTool(~,~,~)%Action with the secundary button
zoom off ; %Desactiva el zoom
[ThHVid,ExHVid,ThHVBM,auxMDLid]=GetV('ThHVid','ExHVid','ThHVBM','auxMDLid');
auxMDLid.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));
ThHVid.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));
ExHVid.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));
ThHVBM.Parent.ButtonDownFcn=@(hObject,eventdata)FHV('ClickAction',hObject,eventdata,guidata(hObject));

function HandleMouseMenu(source,~)
global keep_hv ground_par keep_poi LockAxis change1thickness dataFHV pol_backup
[h,spMDL,auxMDLid,ShowAuxMDL,bmodel]=GetV('h','spMDL','auxMDLid','ShowAuxMDL','bmodel');
LockAxisBack=LockAxis;
% % Formato para pintar back model
[BTA,ALFA,DEN,E]=MODELO(bmodel);
switch get(source,'Label')
    case 'Vs Profile' % Caso Vs
        [source.Parent.Children(2).Checked,source.Parent.Children(4).Checked]=deal('off');
        [source.Parent.Children(3).Checked,spMDL.XLabel.String,spMDL.Parent.Title,ground_par]...
            =deal('on','Velocity Vs [m/s]','Vs Profile',3);
        if ShowAuxMDL % Existencia del back model
            [auxMDLid.XData,auxMDLid.YData(1:end-1)]=deal(BTA,E(1:end-1));
        end

    case 'Vp Profile' % Caso Vp
        [source.Parent.Children(4).Checked,spMDL.XLabel.String,spMDL.Parent.Title,ground_par]...
            =deal('on','Velocity Vp [m/s]','Vp Profile',2);% Asignando variables
        [source.Parent.Children(2).Checked,source.Parent.Children(3).Checked]=deal('off');% Apariencia de la GUI
        if ShowAuxMDL% Existencia del back model
            [auxMDLid.XData,auxMDLid.YData(1:end-1)]=deal(ALFA,E(1:end-1));
        end
        axis auto;
    case 'Density Profile' % Caso de la densidad
        [source.Parent.Children(2).Checked,spMDL.XLabel.String,spMDL.Parent.Title,ground_par]...
            =deal('on','Density [kg/m3]','Density Profile',4);
        [source.Parent.Children(3).Checked,source.Parent.Children(4).Checked]=deal('off');
        if ShowAuxMDL% Existencia del back model
            [auxMDLid.XData,auxMDLid.YData(1:end-1)]=deal(DEN,E(1:end-1));
        end
        axis auto;
    case 'Change one thickness' % Caso del espesor
        [source.Parent.Children(6).Checked,source.Parent.Children(5).Checked,change1thickness]...
            =deal('on','off',1);
    case 'Move one interface'% Caso de la interface
        [source.Parent.Children(5).Checked,source.Parent.Children(6).Checked,change1thickness]...
            =deal('on','off',0);
    case 'Keep Poisson ratio' % Caso Del coediciente de Poisson
        if strcmp(source.Parent.Children(7).Checked,'on')
            [source.Parent.Children(7).Checked,keep_poi]=deal('off',false);
        else
            [source.Parent.Children(7).Checked,keep_poi]=deal('on',true);
        end
    case 'Lock X Axis'
        if LockAxisBack==1,
            [source.Parent.Children(9).Checked,LockAxisBack]=deal('off',0);
        elseif LockAxisBack==2,
            [source.Parent.Children(9).Checked,source.Parent.Children(8).Checked, LockAxisBack]...
                =deal('on','off',1);
        else
            [source.Parent.Children(9).Checked, LockAxisBack]=deal('on',1);
        end
    case 'Lock Y Axis'
        if LockAxisBack==2,
            [source.Parent.Children(8).Checked, LockAxisBack]=deal('off',0);
        elseif LockAxisBack==1,
            [source.Parent.Children(8).Checked,source.Parent.Children(9).Checked,LockAxisBack]...
                =deal('on','off',2);
        else
            [source.Parent.Children(8).Checked, LockAxisBack]=deal('on',2);
        end
    case 'Keep HV ratio'
        if strcmp(source.Parent.Children(1).Checked,'off')
            source.Parent.Children(1).Checked='on';
            keep_hv=true;
        else
            source.Parent.Children(1).Checked='off';
            keep_hv=false;
        end
end
% % Formato para pintar back model
[B,A,R,E]=MODELO(dataFHV(:,1:4));
% Actualizando grafica del perfil
E(end)=sum(dataFHV(:,1))*1.25;
ProfilePolys=[E;A;B;R].';
pol_backup =[ProfilePolys(:,ground_par) ProfilePolys(:,1) ];
setConstrainedPosition(h,pol_backup);

% Actualizando limites de la gráfica del perfil
[bmodel,limdd]=GetV('bmodel','limdd');
if length(bmodel)==1,% bmodel no existe
    source.Parent.Parent.CurrentAxes.XLim=[min(dataFHV(:,ground_par))-min(dataFHV(:,ground_par))*0.05...
        max(dataFHV(:,ground_par))+max(dataFHV(:,ground_par))*0.05  ];
    source.Parent.Parent.CurrentAxes.YLim=[0 ,max(limdd*1.25, E(end)) ];
else
    source.Parent.Parent.CurrentAxes.XLim=[min(min(bmodel(:,ground_par))-min(bmodel(:,ground_par))*0.05 , min(dataFHV(:,ground_par))-min(dataFHV(:,ground_par))*0.05)...
        max(max(bmodel(:,ground_par))+max(bmodel(:,ground_par))*0.05 , max(dataFHV(:,ground_par))+max(dataFHV(:,ground_par))*0.05)];
    source.Parent.Parent.CurrentAxes.YLim=[0 ,max(limdd*1.25, E(end)) ];
end
source.Parent.Parent.CurrentAxes.Children(end).YData(length(E))= max(limdd*1.25, E(end));
LoadV('ZoomOn',0);%Zoom desactivado
LockAxis=LockAxisBack;

%% Import model from *.txt
function LM_Callback(~, ~, handles)
global dataFHV
%Leyendo archivo
[FileName, Path]=uigetfile({'*.txt'},'Load model');
if FileName~=0
    fid=fopen([Path FileName]);
    datos=cell2mat(textscan(fid,' %f %f %f %f ','commentstyle','#'));
    fclose(fid);
    NCF=datos(1,1);
    if isnan(NCF)
        errordlg('Input must be a number','Error');
        return
    elseif NCF<1 || not(mod(NCF,1))==0;
        errordlg('The number of layers should be greater or equal to 1 and integer','Error');
        return
    else
        % Ordenando valores del archivo
        datos(:,5)=(0.5-((datos(:,3))./(datos(:,2))).^2)./(1-((datos(:,3))./(datos(:,2))).^2);
        datos(end,1)=0;%Half space
        [handles.edit1.String,handles.uitable1.Data]=deal(NCF,num2cell(datos(2:end,:)));
        dataFHV=datos(2:end,:);
        % Apariencia de la GUI
        [handles.uitable1.Enable,handles.pushbutton1.Enable,handles.fmin.Enable,handles.fmax.Enable,...
            handles.samples.Enable,handles.bwi.Enable,handles.L.Enable,handles.R.Enable,handles.SM.Enable,...
            handles.LMHV.Enable,handles.uipanel14.Children.Visible,handles.uipanel12.Children(:).Enable]=deal('On');
        handles.LMHV.Value=0;
        LoadV('bmodel',datos(2:end,1:4),'NCFHV',NCF,'ThHVBMLegend',FileName);% Cargando valores al workspace
    end
    auxMDLid=GetV('auxMDLid');
    if length(auxMDLid.XData)~=1,% El molelo auxiliar estaba inicializado, destruyo poligono
        g=GetV('g');
        delete(g.Parent.Children(1));
        MakeLegendHV;
    end
    %LLamando a rutina para perfil interactivo
    VisualFit_v8(auxMDLid,handles);
    g=evalin('base','g');
    % borra el menu del botón derecho
    for ii=1:2*NCF
        delete(g.Children(ii).UIContextMenu.Children)
    end
    % Propiedades de linea para el perfil del backmodel
    delete(g.Children(end).UIContextMenu.Children)
    g.Parent.Children(end).LineStyle='-';
    g.Parent.Children(end).Color=[0 0 0];
    LoadV('limdd',sum(datos(2:end,1)),'ShowAuxMDL',1);% Cargando variables
end

%% H/V for Loaded model
function LMHV_Callback(~, ~,handles) %
[handles.LMHV.Enable,handles.uipanel11.Children(end).Visible]=deal('on');
ThHVBM=GetV('ThHVBM');
if handles.LMHV.Value
    [handles.pushbutton1.Enable,handles.save.Enable]=deal('off');
    CalHV(handles,3)
    handles.pushbutton1.Enable='on';
else
    [ThHVBM.XData,ThHVBM.YData]=deal(nan);
    MakeLegendHV;
end

%% Export model to *.txt
function SM_Callback(~, ~,~ )
global dataFHV
[file,path] = uiputfile('Model.txt','Save Model');
if file~=0
    tex=[path,file];
    fileID = fopen(tex,'w');
    fprintf(fileID,'%s\r' ,num2str(size(dataFHV,1)));
    for row = 1:size(dataFHV,1)
        fprintf(fileID,'%s\r' ,num2str(dataFHV(row,1:4)));
    end
    fclose(fileID);
end
%% Open HV data curve
function LHVD_Callback(~, ~, handles)
[HVFile, Path]=uigetfile({'*.txt'},'Open HV Curve (Data)');
if HVFile~=0
    fid=fopen([Path HVFile]);
    HVD=(textscan(fid,' %f %f %f ','delimiter','\t','commentstyle','#'));% Read file  delimited simple tabulador
    fclose(fid);
    if length(HVD{1})<=1
        fid=fopen([Path HVFile]);
        HVD=(textscan(fid,' %f %f %f ','delimiter',',','commentstyle','#'));% Read file delimited with simple comma
        fclose(fid);
    end
    if length(HVD{1})<=1
        fid=fopen([Path HVFile]);
        HVD=(textscan(fid,' %f %f %f ','delimiter',' ','commentstyle','#'));% Read file delimited with simple space
        fclose(fid);
    end
    if isempty(HVD{1}) || length(HVD{1})<=1  % The file is empty or headers uncomment
        warndlg('the file is empty or the headers are not comment with ´´#´. Check the file ','!! Warning !!')
        pass=false;% Cancela la operación
    else
        if  ~isnan (HVD{3}(1,1)) % Si la std no esta vacia
            if length(HVD{1})==length(HVD{2}) && length(HVD{1})==length(HVD{3}) % Vectores del mismo tamaño
                pass=true;
            elseif isnan (HVD{3}(1,1)) % cuando std esta vacio
                if  length(HVD{1})==length(HVD{2})
                    pass=true;
                else
                    pass=false;
                    warndlg('the file is empty or the headers are not comment with ´´#´. Check the file ','!! Warning !!')
                end
            else
                pass=false;
                warndlg('The data are not match. Check the file ','!! Warning !!')
            end
        else
            if length(HVD{1})==length(HVD{2})
                pass=true;
            else
                pass=false;
                warndlg('The data are not match. Check the file ','!! Warning !!')
            end
        end
    end
    if pass
        ClearHVDataOption=['Clear ' HVFile];
        [ExHVid,cF]=GetV('ExHVid','cF');
        if length(cF.Children)==5, % Borrar del menu el rotulo de borrado de un posible antiguo modelo
            delete(cF.Children(1));
        end
        uimenu(cF,'Label',ClearHVDataOption,'Callback',@setlinestyle);
        [ExHVid.XData,ExHVid.YData,ExHVid.LData,ExHVid.UData]=deal(nan);
        [ExHVid.XData,ExHVid.YData]=deal(HVD{1},HVD{2});
        if length (HVD{3}) == length (HVD{2})
            [ExHVid.LData,ExHVid.UData]=deal(HVD{3});
        end
        LoadV('ExHVLegend',HVFile,'ClearHVDataOption',ClearHVDataOption,'HVFHVD',HVD{2},'FFHVD',HVD{1})% Load variables workspace
        MakeLegendHV;% Leyenda al eje de los HV
        handles.uipanel11.Children(end).Visible='on';%GUI Changes appearance
    end
end

%% Calculate HV Curve from HVf.exe
function CalHV(handles,idHV)
global dataFHV
pas=1;% id control
if idHV==1 % Id theorical HV
    id='ThHVid';
else
    id='ThHVBM'; % Id back model
    data =GetV('bmodel'); % back model
end
%Get variables
[HVid,apsvFHV,NC,fmin,fmax,R,L,nksBW,n,IsLogF]=...
    GetV(id,'apsvFHV','NCFHV', 'fminFHV','fmaxFHV','RFHV', 'LFHV','nksBW','nFHV','IsLogF');
[handles.save.Enable,handles.pushbutton1.Enable,handles.uipanel11.Children(end).Visible]=deal('off', 'off','on');
% Condicionantes Ramas surface waves and body
if R==0 && L~=0 && nksBW==0
    errordlg('Impossible to obtain the ratio H/V','Error')
elseif R==0 && L==0 && nksBW==0
    errordlg('Impossible to obtain the ratio H/V','Error')
else
    % Creacion de modelo *.txt
    if idHV==1
        par=([[NC; dataFHV(:,1)]  [0 ;dataFHV(:,2) ] [0 ;dataFHV(:,3) ] [0 ;dataFHV(:,4) ]]); % Modelo Test HV
    else
        par=([[NC;data(:,1)] [0;data(:,2)] [0;data(:,3)] [0;data(:,4)]]);% backmodel HV
    end
    % Function of execute HVf.exe
    [HVFHV,FFHV]=HVPD(par,fmin,fmax,n,R,L,nksBW,IsLogF,apsvFHV);
    % Condicionantes
    if isnan(sum(HVFHV)) || isinf(sum(HVFHV))  ;
        imal=isnan(HVFHV) | isinf(HVFHV);
        if length(imal)>length(HVFHV)*0.1
            errordlg('Impossible to obtain the H/V ratio','Error')
            pas=0;% id control
        else
            HVFHV(imal)=interp1(A.Fobs(~imal),HVFHV(~imal),A.Fobs(imal),'linear','spline'); % Interpolación
            pas=1; % id control
        end
    end
    % GUI Graphic to curve
    if pas==1% id control
        [HVid.XData,HVid.YData]=deal(FFHV,HVFHV);% Grafica en blanco
        MakeLegendHV;% Update Legend
        if idHV==1% solo guarda el HV test
            LoadV('HVFHV',HVFHV,'FFHV',FFHV); % Save data HV (frequency, Amplitude)
        end
    end
end
[handles.save.Enable,handles.pushbutton1.Enable]=deal('on');

%% Function for Legend for HV axes
function MakeLegendHV()
[spHV,ExHVLegend,ThHVBMLegend]=GetV('spHV','ExHVLegend','ThHVBMLegend');
valid_oids=[];
valid_labels={};
if isempty(find(isnan(spHV.Children(1).XData),1))
    [valid_oids,valid_labels{end+1}]=deal([valid_oids,spHV.Children(1)],ThHVBMLegend);
end
if isempty(find(isnan(spHV.Children(2).XData),1))
    [valid_oids,valid_labels{end+1}]=deal([valid_oids,spHV.Children(2)],ExHVLegend);
end
if isempty(find(isnan(spHV.Children(3).XData),1))
    [valid_oids,valid_labels{end+1}]=deal([valid_oids,spHV.Children(3)],'H/V test model');
end
if ~isempty(valid_oids)
    legend(spHV,valid_oids,valid_labels);
    legend(spHV,'show');
else
    legend(spHV,'hide');
end

%% GUI Function for Visualization of HV subplot
function setlinestyle(source,~)
[spHV,ClearHVDataOption,cF]=GetV('spHV','ClearHVDataOption','cF');
switch get(source,'Label')
    case 'XScale: Log'
        spHV.XScale='log';
    case 'XScale: Linear'
        spHV.XScale='linear';
    case 'YScale: Log'
        spHV.YScale='log';
    case 'YScale: Linear'
        spHV.YScale='linear';
    case ClearHVDataOption
        [ExHVid]=GetV('ExHVid');
        [ExHVid.XData,ExHVid.YData,ExHVid.LData,ExHVid.UData]=deal(nan);
        delete(cF.Children(1));
        MakeLegendHV;
        %         if ~isempty(find(~isnan(FFHV),1)),xlim(spHV,[min(FFHV) max(FFHV)]);end
        %         if ~isempty(find(~isnan(HVFHV),1)),ylim(spHV,[min(HVFHV)-min(HVFHV)*0.5 max(HVFHV)+max(HVFHV)*0.1]);end
        LoadV('HVFHVD',nan,'FFHVD',nan,'ExHVLegend','');
end
%% GUI close interface option
function figure1_CloseRequestFcn(~, ~, ~)
choice = questdlg('Quitting?','Are you Sure?','Yes','No','No');
switch choice
    case 'Yes'
        closereq
    otherwise
        return
end

function Exit_Callback(~, ~, ~)
choice = questdlg('Quitting?','Are you Sure?','Yes','No','No');
switch choice
    case 'Yes'
        closereq
    otherwise
        return
end

function Menu_Callback(~, ~, ~)
function edit1_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function varargout = FHV_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
function fmin_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function fmax_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function samples_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function bwi_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function L_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function R_CreateFcn(hObject, ~, ~)
Auxiliary(hObject)
function Auxiliary(hObject)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

