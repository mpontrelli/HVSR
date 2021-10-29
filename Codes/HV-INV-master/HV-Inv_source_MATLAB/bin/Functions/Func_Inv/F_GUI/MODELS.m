%     Copyright (C) 2014,2016 José Piña-Flores, Antonio García-Jerez.
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License version 3 as
%     published by the Free Software Foundation.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Funcion GUI menu MODELS
% Generating H/V and Dispersion graphics curves, result the inversion

%% Function for MODELS GUI menu
function varargout = MODELS(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @MODELS_OpeningFcn, ...
    'gui_OutputFcn',  @MODELS_OutputFcn, ...
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

%% Module Input the information by user
function MODELS_OpeningFcn(hObject, ~, handles, varargin)
L2=GetV('L2sens'); % Norm L2
[a2,  ~]=sort(L2,'descend');
bm=ceil(a2(1)*1000)/1000;
handles.erri.String=['The minimum misfit is',' ',num2str(a2(end))];
handles.edit1.String=num2str(bm);
LoadV('MIS',bm,'ER',a2,'sel',1);
handles.output = hObject;
guidata(hObject, handles);
movegui(handles.figure1,'center')

%% Input information of misfit or number of model for grasphics
function edit1_Callback(hObject, ~, ~)
[ERROR,sel]=GetV('ER','sel');
E=ERROR(end);
MIS = str2double(hObject.String);
if sel==1;
    if MIS<E
        errordlg(['The number of Maximun Misfit must be greater than ', num2str(E)],'Error');
    elseif isnan(MIS)
        errordlg('Input must be a number','Error');
    end
else
    if MIS<1
        errordlg(['The number of Models must be greater than ', num2str(1)],'Error');
    elseif isnan(MIS)
        errordlg('Input must be a number','Error');
    end
end
LoadV('MIS',MIS);

function pushbutton1_Callback(~, ~, ~)
close
[sel,MIS,flag,a,L2,NC,stdgHV,stdgDC]=GetV('sel','MIS','flag','asens','L2sens','NC','stdgHV','stdgDC');
[a2, b2]=sort(L2,'descend');
if sel==1;
    pro=find((MIS>a2)==1);
else
    if MIS>=length(a2)
        pro=find((a2(1)>=a2)==1);
    else
        pro1=a2(end-(MIS-1):end);
        pro=find((pro1(1)>=a2)==1);
    end
end
%% Elastic parameters
dg=length(pro);
[E,ALFA,BTA,RO]=deal(zeros(dg,NC*2));
% [L2bar,modeltex]=deal(zeros(dg,1),zeros(NC,4,dg));
% fid=fopen([pwd '\models\models.mout'],'w');
% fprintf(fid,'%s',[' #Thickness ' ' #Vp ' ' #Vs '  ' #Density ']);
% fprintf(fid,'\n');
% mms=waitbar(0,{'Creating Models','Please wait'});
% mms.CloseRequestFcn=' ';
mms = msgbox({'Creating Figures','Please wait'});
for i=1:dg
    [BTA(i,:),ALFA(i,:),RO(i,:),E(i,:)]=MODELO(a(:,:,b2(pro(i))));
    %     [modeltex(:,:,i),L2bar(i)]=deal(a(:,:,b2(pro(i))),L2(b2(pro(i))));
    %waitbar(i/dg);
end
L2bar=L2(b2(pro));
%;

%% Plot velocity profile
fig5=figure(5);
%% Profile Vs
%movegui(fig5,'center');
set(0,'DefaultAxesColorOrder',parula(dg))
% Plot Velocity profile VS
P1=subplot(1,3,2);
P1.NextPlot='add';
plot (P1,BTA',E','LineWidth',3.5);
h2=plot (P1,BTA(end,:),E(end,:),'r','LineWidth',2.5);
[h2.Parent.YDir,h2.Parent.XLabel.String,h2.Parent.YLabel.String]=...
    deal('reverse','$Velocity$ $V_s$ [${m} \over {s}$]','$Depth$ [$m$]');
l2=legend(h2,'Best model V_s','Location','best');
[fig5.CurrentAxes.FontSize,fig5.CurrentAxes.FontUnits]=deal(10,'normalized');
%% Profile Vp
P2=subplot(1,3,1);
P2.NextPlot='add';
plot(P2,ALFA',E','LineWidth',3.5);
% Plot Velocity profile Vp
h1=plot (P2,ALFA(end,:),E(end,:),':r','LineWidth',2.5);
% Settings Plot
[h1.Parent.YDir,h1.Parent.XLabel.String,h1.Parent.YLabel.String]=...
    deal('reverse','$Velocity$ $V_p$ [${m} \over {s}$]','$Depth$ [$m$]');
l1=legend(h1,'Best model V_p','Location','best');
[fig5.CurrentAxes.FontSize,fig5.CurrentAxes.FontUnits]=deal(10,'normalized');
%% Profile Density
P3=subplot(1,3,3);
P3.NextPlot='add';
plot(P3,RO',E','LineWidth',3.5);
% Plot Velocity profile density
h3=plot (P3,RO(end,:),E(end,:),'-r','LineWidth',2.5);
% Settings Plot
[h3.Parent.YDir,h3.Parent.XLabel.String,h3.Parent.YLabel.String]=...
    deal('reverse','$Density$ [${Kg} \over {m^3}$]','$Depth$ [$m$]');
l3=legend(h3,'Best model Rho','Location','best');
[fig5.CurrentAxes.FontSize,fig5.CurrentAxes.FontUnits]=deal(10,'normalized');
[l1.FontSize,l2.FontSize,l3.FontSize]=deal(12);
[h3.Parent.XLabel.Interpreter,h3.Parent.YLabel.Interpreter,h2.Parent.XLabel.Interpreter...
    ,h2.Parent.YLabel.Interpreter,h1.Parent.XLabel.Interpreter,h1.Parent.YLabel.Interpreter]=deal('latex');
% Colorbar
colormap(parula(dg));
if strcmp(version ('-release'),'2015a') || strcmp(version ('-release'),'2015b');
    h=colorbar('Ticks',(linspace(0,1,11)),'YTickLabel',round(L2bar(ceil(linspace(1,dg,11))),2),'Fontsize',12);
else
    h=colorbar('YTick',ceil(linspace(1,dg,11)),'YTickLabel',(L2bar(ceil(linspace(1,dg,11)))),'Fontsize',12);
end
[h.YLabel.String,h.YLabel.FontSize,fig5.Name,fig5.NumberTitle]=deal('Misfit',10,'Models','off');
h.YLabel.FontUnits='normalized';
delete (mms)
%% Plot Curves HV and Dispersion

[C,DCALL,HVALL]=GetV('C','DCALL','HVALL');
%% Plot HV
if flag==1 || flag==3
    mms = msgbox({'Creating Figures H/V','Please wait'});
    HVG=HVALL(:,b2(pro));
    %% Plot HV curves
    fig6=figure(6);
    %movegui(fig6,'center');
    set(0,'DefaultAxesColorOrder',parula(dg))
    semilogx(C.HVFobs,HVG,'LineWidth',5)
    hold on
    h2v=semilogx(C.HVFobs,HVG(:,end),'-r','LineWidth',2.5);
    if stdgHV~=3;
        h1v=errorbar(C.HVFobs,C.HVAobs,C.HVDS,'-k','LineWidth',2.5);
        set(get(h1v,'Children'),{'LineWidth'},{2.5; 0.5})
    else
        h1v=semilogx(C.HVFobs,C.HVAobs','-k','LineWidth',2.5);
    end
    % Settings plot
    [h2v.Parent.XLabel.String,h2v.Parent.YLabel.String, h2v.Parent.XLim,h2v.Parent.Title.String] =...
        deal('$Frequency$ [$Hz$]','$Amplitude$ [${H} \over {V}$]',[min(C.HVFobs) max(C.HVFobs)],'$HV$ $ratio$');
    l4=legend([h2v h1v],'Best HV',C.HVFile{1},'Location','best');
    [l4.Interpreter,l4.FontSize]=deal('none',10);
    [fig6.CurrentAxes.FontSize,fig6.CurrentAxes.FontUnits]=deal(10,'normalized');
    [h2v.Parent.XLabel.Interpreter,h2v.Parent.YLabel.Interpreter,h2v.Parent.Title.Interpreter]=deal('latex');
    % Colorbar
    colormap(parula(dg));
    if strcmp(version ('-release'),'2015a') || strcmp(version ('-release'),'2015b');
        hv=colorbar('Ticks',(linspace(0,1,11)),'YTickLabel',round(L2bar(ceil(linspace(1,dg,11))),2),'Fontsize',12);
    else
        hv=colorbar('YTick',ceil(linspace(1,dg,11)),'YTickLabel',(L2bar(ceil(linspace(1,dg,11)))),'Fontsize',12);
    end
    [hv.YLabel.String,hv.YLabel.FontSize,fig7.Name,fig7.NumberTitle,hv.YLabel.Units]=...
        deal('Misfit',10,'DC curves','off','normalized');
    delete(mms)
end

%% Plot Dispersion Curve
if flag==2 || flag==3
    [vel,pol]=GetV('VEL','POL');
    mms = msgbox({'Creating Figures DC','Please wait'});
    CG=DCALL(:,b2(pro));
    %% Plot dispersion curves
    fig7=figure(7);
    %movegui(fig7,'center');
    set(0,'DefaultAxesColorOrder',parula(dg))
    semilogx(C.DCFobs,CG,'LineWidth',5)
    hold on
    h2c=semilogx(C.DCFobs,CG(:,end),'-r','LineWidth',2.5);
    if  stdgDC~=3;
        h1c=errorbar(C.DCFobs,C.DCVobs,C.DCDS,'-k','LineWidth',2.5);
        set(get(h1c,'Children'),{'LineWidth'},{2.5; 0.5})
    else
        h1c=semilogx(C.DCFobs,C.DCVobs,'-k','LineWidth',2.5);
    end
    % Settings plot
    if vel==1
        h2c.Parent.YLabel.String= '$Phase$ $velocity$ [${m} \over {s}$]';
    else
        h2c.Parent.YLabel.String ='$Group$ $velocity$ [${m} \over {s}$]';
    end
    
    if pol==1
    h2c.Parent.Title.String='$Rayleigh$ $waves$';
    else
    h2c.Parent.Title.String='$Love$ $waves$';
    end
    [h2c.Parent.XLabel.String, h2c.Parent.XLim] =deal('$Frequency$ [$Hz$]',[min(C.DCFobs) max(C.DCFobs)]);
    l4=legend([h1c h2c], C.DCFile{1},'Best DC','Location','best');
    [l4.Interpreter,l4.FontSize]=deal('none',10);
    [fig7.CurrentAxes.FontSize,fig7.CurrentAxes.FontUnits]=deal(10,'normalized');
    [h2c.Parent.XLabel.Interpreter,h2c.Parent.YLabel.Interpreter,h2c.Parent.Title.Interpreter]=deal('latex');
    colormap(parula(dg));
    if strcmp(version ('-release'),'2015a') || strcmp(version ('-release'),'2015b');
        hc=colorbar('Ticks',(linspace(0,1,11)),'YTickLabel',round(L2bar(ceil(linspace(1,dg,11))),2),'Fontsize',12);
    else
        hc=colorbar('YTick',ceil(linspace(1,dg,11)),'YTickLabel',(L2bar(ceil(linspace(1,dg,11)))),'Fontsize',12);
    end
    [hc.YLabel.String,hc.YLabel.FontSize,fig7.Name,fig7.NumberTitle,hc.YLabel.Units]=...
        deal('Misfit',10,'DC curves','off','normalized');
    delete(mms)
end
hg=[h1 h2];
LoadV('hg',hg);
modelo_medio4;


%% Function for closed gui
function pushbutton2_Callback(~, ~, ~)
close

%% Funtion for change the option for misfit or number of models
function uibuttongroup1_SelectionChangedFcn(hObject,~, handles)
L2=GetV('L2sens');
[a2,  ~]=sort(L2,'descend');
%% Option for misfit value
if  hObject==handles.uno
    sel=1;
    bm=ceil(a2(1)*1000)/1000;
    handles.erri.String=['The minimum misfit is',' ',num2str(a2(end))];
    handles.edit1.String=num2str(bm);
    %% option for number of models
elseif hObject==handles.dos
    sel=2;
    bm=length(a2);
    handles.erri.String=['The maximum number of models is',' ',num2str(bm)];
    handles.edit1.String=num2str(bm);
end
LoadV('MIS',bm,'ER',a2,'sel',sel);

function erri_CreateFcn(~, ~, ~)
function edit1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function varargout = MODELS_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
%% Lines for write file whit models information
% fid=fopen([pwd '\models\HVout.out'],'w');
%     fprintf(fid,'%s',[' #Frequency ' ' Amplitude HV ']);
%     fprintf(fid,'\n');
%     fprintf(fid,'%s \b','#Misfit  ');
%     fprintf(fid,'%10.5f \b',L2bar');
%     fprintf(fid,'\n');