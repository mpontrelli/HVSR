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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Function  for input for initial model by user
%% Function for initialmodel GUI

function varargout = initialmodel(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @initialmodel_OpeningFcn, ...
    'gui_OutputFcn',  @initialmodel_OutputFcn, ...
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

%% GUI Function elastic parameters by default
function initialmodel_OpeningFcn(hObject, ~, handles, varargin)
% GUI Variables
handles.output = hObject;
guidata(hObject, handles);
movegui(handles.figure1,'center');
%% Elastic parameters
NC=GetV('NC');
if NC~=0
    [data]=GetV('data');
    ini1=data(:,[1 3 5 7])+(-data(:,[1 3 5 7])+data(:,[2 4 6 8]))./2;
    ini1(end,2:3)=ini1(end,2:3)*1.1;
    [ini1(:,5),nu1F]=deal((0.5-(ini1(:,3)./ini1(:,2)).^2)./(1-(ini1(:,3)./ini1(:,2)).^2));
    %Edit table of model parameters
    handles.uitable2.Data=num2cell(ini1);
    LoadV('ini1d',ini1,'nu1inv',nu1F,'CHT',0);
else
    [handles.uitable2.Enable,handles.pushbutton1.Enable]=deal('off');
    handles.text2.String='Input the Initial model from file *.txt';
end

%% GUI Function accept initial model by default
function pushbutton1_Callback(~, ~,~)
[CHT,ini1,handles2]=GetV('CHT','ini1d','handles');
if CHT
    %% GUI
    [handles2.START.Visible,handles2.TableParameters.Enable,...
        handles2.SPara.Enable,handles2.LV.Enable,handles2.T0.Enable,...
        handles2.Redu.Enable,handles2.Inversiont.Enable]=deal('On');
    %% Set variables
    NC=size(ini1,1);
    handles2.edit1.String=NC;
    data=zeros(NC,10);
    [data(:,1),data(:,2),data(:,3),data(:,4),data(:,5),data(:,6),data(:,9),data(:,10)]=...
        deal(ini1(:,1)-ini1(:,1).*0.5,ini1(:,1)+ini1(:,1).*0.5,ini1(:,2),ini1(:,2),ini1(:,3),ini1(:,3),0.1,0.49);
    [data(end,1),data(end,2)]=deal(0);
    [data(:,7),data(:,8)]=deal(ini1(:,4));% min and max density
    handles2.TableParameters.Data=num2cell(data);
    var=zeros(NC,15);
    [var(:,1:3:13),var(:,2:3:14)]=deal(data(:,1:2:9),data(:,2:2:10));% Minima to matrix var% Maxima to matrix var
    var(:,3:3:15)=var(:,2:3:14)-var(:,1:3:13);
    limdd=sum(var(:,2)); %max. depth down to the halfspace
    LoadV('var',var,'limdd',limdd,'data',data,'NC',NC,'ini1',ini1(:,1:4));
else
    LoadV('ini1',ini1(:,1:4));
end
close;

%% GUI Function - Cancel initial model
function pushbutton2_Callback(~, ~,~)
close;

%% GUI Function edit elastic parameters by user (Table parameters)
function uitable2_CellEditCallback(hObject,  ~, handles)
dataF = cell2mat (hObject.Data);
nu1F=GetV('nu1inv');
%% Verifying elastic parameters depending of range
if sum(find(nu1F~=dataF(:,5)))>=1
    dataF(:,2)=real(dataF(:,3).*(sqrt((2*(1-dataF(:,5)))./(1-(2.*dataF(:,5))))));
else
    dataF(:,5)=(0.5-(dataF(:,3)./dataF(:,2)).^2)./(1-(dataF(:,3)./dataF(:,2)).^2);
end
dataF(end,1)=0;
handles.uitable2.Data=num2cell(dataF);
LoadV('ini1d',dataF,'ini1',dataF(:,1:4),'nu1inv',dataF(:,5),'CHT',1);


function pushbutton3_Callback(~, ~,handles)
[FileName, Path]=uigetfile({'*.txt'},'Open Parameters');
if FileName~=0
    fid=fopen([Path FileName]);
    datos=textscan(fid,' %f %f %f %f %f  ','commentstyle','#');
    fclose(fid);
    NC=datos{1}(1);
    ini1d=cell2mat(datos);
    ini1=ini1d(2:end,:);
    %% Verifying elastic parameters
    ini1(:,5)=(0.5-(ini1(:,3)./ini1(:,2)).^2)./(1-(ini1(:,3)./ini1(:,2)).^2);
    ini1(end,1)=0;
    LoadV('nu1inv',ini1(:,5),'ini1d',ini1);
    [handles.uitable2.Enable,handles.pushbutton1.Enable]=deal('on');
    CHT=0;
    [NCT]=GetV('NC');% Numero de capas de la tabla de intervalos  (NC=0 no existe)
    if NCT~=0
        if NCT==NC
            handles.uitable2.Data=num2cell(ini1);
        else
            choice = questdlg('The number of layers doesn`t match the model parameters table.The table will be restarted.','Warning','Ok','Cancel','Ok');
            % Handle response
            switch choice
                case 'Ok'
                    CHT=1;% Autorizacion de actualizacion de la tabla de intervalos CHT=1;
                    handles.uitable2.Data=num2cell(ini1); %% Escribe en la tabla pequeña
                    LoadV('numele1',ini1);
                case 'Cancel'
                    CHT=0;% Negación de actualizacion de la tabla de intervalos CHT=0;
            end
        end
    else
        CHT=1;% Autorizacion de actualizacion de la tabla de intervalos CHT=1;
        handles.uitable2.Data=num2cell(ini1); %% Escribe en la tabla pequeña
        LoadV('numele1',ini1);
    end
    LoadV('CHT',CHT)
    handles.uitable2.ColumnEditable=true(1,5);
end

%% Additional GUI Functions
function varargout = initialmodel_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
