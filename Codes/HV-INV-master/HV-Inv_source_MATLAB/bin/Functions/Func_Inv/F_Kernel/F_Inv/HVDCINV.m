%     Copyright (C) 2014,2016 José Piña-Flores, Antonio García-Jerez.
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License version 3 as
%     published by the Free Software Foundation..
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%      #ok<*DEFNU>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for inversion methods and GUI control

function HVDCINV(handles,ER,Cont)
%INPUT
%handles elements of GUI
%ER struc for information during inversion
%Cont contador de iteraciones
%%  Getting variables from the workspace
% A%Information of iterations
% A.open% Enable-disable parallel
% G.Pl%Tyoe axes control
% C=Information of targets curves
% A.HHS Hard Half Space option (1=yes or 0=no)
% A.apsv Control of imaginary frequency to stabilize body-wave integrals
% pmy=% pmy == 1 to equalize by the number of samples
% C.ww= % weight of dispersion curve in joint inversion
% ini1=% initial model (for some methods only)
% A.HHS=;%Halfspace
% A.ZLVs=;%Zone Low Velocity
% A.ZLVp=%Zone Low Velocity
% A.NC=%Number of layers
% A.var=%Range of parameters
% A.T0=%Temperature
% A.EA=;%Error aceptation (SA MSA)
% A.AP=%aceptaction Probability  (SA MSA)
% A.RT=%Reduce temperature
% Tinv=%Type of inversion (HV DC or HVDC)
%Get Variables from workspace
[A,A.open,G.Pl,C,A.HHS,A.apsv,pmy,C.ww,ini1,A.HHS,A.ZLVs,A.ZLVp,A.NC,A.var,A.T0,A.EA,A.AP,...
    A.RT,Tinv]=GetV('A','open','Pl','C','HHS','apsv','pmy','ww','ini1','HHS','ZLVs','ZLVp','NC','var', 'T0','EA','AP', 'RT','Tinv');
%% Calculation of weights for joint inversion when equalization is required
if pmy &&  ~isempty(C.HVFobs) && ~isempty(C.DCFobs)
    % The weight for the dispersion curve will be:
    C.ww=C.nmHV/(C.nmDC+C.nmHV);% = const/C.nmDC, with const=C.nmDC*C.nmHV/(C.nmDC+C.nmHV)
    % The weight for H/V will be (1-C.ww) = const/C.nmHV
end
%% flag for variables HV curve
if  ~isempty(C.HVFobs)
    %Informacion de la diferente contribucion de ondas para el calculo del
    %HV
    [C.NR,C.NL,C.dk,C.mdk]=GetV('ramasR','ramasL','kb_by_dk','mkb_by_dk');
end
%% flag for variables DC curve
if  ~isempty(C.DCFobs)
    % Obteniendo informacion del tipo de curva de dispersión
    [C.POL,C.VEL,C.MODE]=GetV('POL','VEL','MODE');
end
A.ST=A.Nsubiter;% Numero de iteraciones
A.limd= sum(A.var(:,2))+0.1*sum(A.var(:,2));% Limit halfspace
if isempty(ini1)% Si no existe modelo inicial
    INPO=GetV('INPO');% Numero de poblacion inicial
else
    INPO=1;
end
%% Control processes in parallel.
if A.open==0
    MM=1;
else
    if Tinv==3
        MM=GetV('NUMW');% NUmero de nucleos habilitados para procesos en paralelo
    else
        MM=1;% Si no estan habilitados los workers 
    end
end
%% Variable initialization for processes in parallel
if A.open==1
    NPI=floor(INPO/MM);% distribucion del cálculo de las curva dependiendo el numero de workers
else
    NPI=INPO;
end
% HVT=zeros( C.nmHV,MM*NPI);%HV curve
% DCT=zeros( C.nmDC,MM*NPI);%Disp curve
% HV=zeros(C.nmHV,1);%HV curve
% DC=zeros(C.nmDC,1);%Disp curve
% DCTp=zeros(C.nmHV,NPI);%HV curve in parallel
% HVTp=zeros(C.nmDC,NPI);%DC curve in parallel
%modp=zeros(A.NC,4,NPI)*nan; %model in parallel
% L2i=zeros(1,MM*NPI)*nan;% norm L2
% L2ip=zeros(1,NPI)*nan;% norm L2 in parallel
% It=(A.Nsubiter+1)*A.Niter+floor(A.Niterf);% Number of acomulated iterations
% A.flag=zeros(MM,1);% Flag for SAM in parallel
% A.in=zeros(A.NC,4);% initial Model for perturbation
[HVT,DCT,HV,DC,DCTp,HVTp]=deal(zeros( C.nmHV,MM*NPI),zeros( C.nmDC,MM*NPI),...
    zeros(C.nmHV,1),zeros(C.nmDC,1),zeros(C.nmDC,NPI),zeros(C.nmHV,NPI));
[mod,modp]=deal(nan(A.NC,4,MM*NPI),nan(A.NC,4,NPI));
[L2i,L2ip,It,A.flag,A.in]=deal(nan(1,MM*NPI),nan(1,NPI),(A.Nsubiter+1)*A.Niter+floor(A.Niterf),...
    zeros(MM,1),zeros(A.NC,4));
% Initializing auxiliar matrices
% The Vs of the halfspace must be the larger one. This flag can be changed below to false for not implemented cases
if find(A.var(1:A.NC-1,7)>A.var(A.NC,8)),errordlg('Impossible to set hard halfspace condition. Turning that condition to false');uiwait;A.HHS=false;end
if A.ZLVs&&A.ZLVp&&A.HHS % only case with HHS implemented
    [A.pol,A.limseg]=deal([]);
    [~,~,~,A.polHHS,A.facHHS,A.ProbLimSupSegHHS,A.LimVelSegHHS]=ini1_plus(A.NC,A.var,A.ZLVs,A.ZLVp,[],[],A.HHS);
else
    A.HHS=false;[A.polHHS,A.facHHS,A.ProbLimSupSegHHS,A.LimVelSegHHS]=deal([]);
    [~,A.pol,A.limseg]=ini1_plus(A.NC,A.var,A.ZLVs,A.ZLVp);
end
pass=1;
%% Waiting Message
if isempty(ini1);
    mms=waitbar(0,{'Initial random exploration of the parameters space','Please wait'},'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
    setappdata(mms,'canceling',0);
    set(mms, 'CloseRequestFcn','')
end
handles.STOP.Enable='off';
%% Generation of random models (random search)
if A.open==0
    for i=1:INPO;
        if isempty(ini1)
            if getappdata(mms,'canceling') % Stopped by user
                [handles.STOP.UserData, handles.START.Visible]=deal(1,'on');
                delete(mms);
                pass=0;
                break;
            end
            [mod(:,:,i),HV,DC,~,L2i(:,i)]=CTF(A,C,i,1,G);
            [DCT(:,i),HVT(:,i)]=deal(DC,HV);
            waitbar(i/INPO);
        else
            [mod(:,:,i),HV(:,i),DC(:,i),~,L2i(:,i)]=CTF(A,C,i,3,G);
            [DCT(:,i),HVT(:,i)]=deal(DC,HV);
            if isnan(L2i);
                errordlg(sprintf('Forward calculation failed for the initial model\n\nPlease, enter a new initial model'),'Error');
                [handles.STOP.UserData, handles.START.Visible]=deal(1,'on');
                pass=0;
            end
        end
    end
elseif A.open~=0 && isempty(ini1); % Modo de calculó de curvas en paralelo (Initial random)
    kij=0;
    for ii=1:MM
        % cancelación de proceso del usuario
        if getappdata(mms,'canceling')
            [handles.STOP.UserData, handles.START.Visible]=deal(1,'on');
            delete(mms);
            pass=0;
            break;
        end
        % processes in parallel
        parfor i=1:floor(INPO/MM);
            [modp(:,:,i),HV,DC,~,L2ip(:,i)]=CTF(A,C,i,1,G);
            [DCTp(:,i),HVTp(:,i)]=deal(DC,HV);
        end
        mod(:,:,((kij+1):(ii*floor(INPO/MM))))=modp;
        L2i(((kij+1):(ii*floor(INPO/MM))))=L2ip;
        % Cálculo para las curvas de dispersión modelos y misfit
        if  ~isempty(C.DCFobs)
            DCT(:,((kij+1):(ii*floor(INPO/MM))))=DCTp;
        end
        %Calculo para las curvas HV modelos y misfit
        if  ~isempty(C.HVFobs)
            HVT(:,((kij+1):(ii*floor(INPO/MM))))=HVTp;
        end
        kij=ii*floor(INPO/MM);
        waitbar(ii/MM);
    end
elseif A.open~=0 && ~isempty(ini1); %Cuando existe un modelo inicial se calcula el misfit 
    [mod(:,:,1),HV,DC,~,L2i(:,1)]=CTF(A,C,1,3,G);
    [DCT(:,1),HVT(:,1)]=deal(DC,HV);
    if isnan(L2i);
        errordlg(sprintf('Forward calculation failed for the initial model\n\nPlease, enter a new initial model'),'Error');
        [handles.STOP.UserData, handles.START.Visible]=deal(1,'on');
        pass=0;
    end
end
%% Optimal allocation model of Grid-search and parameter settings

if isempty(ini1);
    delete (mms);
end
handles.STOP.Enable= 'on';% Habilitacion del boton STOP
Cont.all=0;
[~ ,b]=sort(L2i);%Ordenando el misfit de menor a mayor
aopt=mod(:,:,b(1));%Ordenando midelos de menor a mayor el misfit
[L2i1,A.L2opt,A.L2optimo]=deal(L2i(b(1)));
[A.NN1,A.NN2]=deal(A.NN/100);
[DCopt,HVopt,A.aini]=deal(DCT(:,b(1)),HVT(:,b(1)),aopt);
% L2i1=L2i(b(1));% Auxiliar norm L2
% A.L2opt=L2i(b(1));% Best norm L2
% A.L2optimo=L2i(b(1));% Best misfit model
% DCopt=DCT(:,b(1));% Best disp curve
% HVopt=HVT(:,b(1));% Best HV curve
% A.aini=aopt;% initial model
% A.NN1=A.NN/100;% Perturbation of range in percentage
% A.NN2=A.NN/100;% Perturbation of range in percentage
if ~isempty(ini1)
    [modini1,HVTini1,DCTini1,L2ini1]=deal(mod(:,:,1),HVT(:,1),DCT(:,1),L2i(:,1));
end
handles.textmisfit.String=['#Evaluated models',' ','1',' ','Min. Misfit',' ',num2str(A.L2optimo)];
%% GUI graphing module
% GUI graphical error
[F1,F2,Perr]=GetV('F1','F2','Perr');
if pass==1
    if isempty(F1);
        handles.ploterror.NextPlot='add';
        [Perr,F1,F2]=plotyy(handles.ploterror,0,[0; 0],0,0);
        uistack(Perr(1),'top');
        %% Graphic settings %apariencia de la grafica del error
        % Escalas automaticas
        [Perr(2).YScale,Perr(1).YScale]=deal('Log');
        [F1(2).Parent.YTickMode,Perr(1).YLimMode,Perr(1).XLimMode,Perr(1).XTickMode,Perr(2).YLimMode,...
            Perr(2).XLimMode,Perr(2).XTickMode,handles.ploterror.XLimMode,handles.ploterror.YLimMode]=deal('auto');
        % Legendas de ejes
        Perr(1).YLabel.String='Misfit';
        if Tinv==2 || Tinv==3
            Perr(2).YLabel.String='Temperature';
        elseif Tinv==5 || Tinv==4 || Tinv==1
            [Perr(2).YLabel.String,Perr(2).YTickMode,Perr(2).YTick]=deal(' ','manual',[]);
        end
        %Legend error graphic
        fgd1=legend(handles.ploterror,'Minimum misfit','Current misfit');
        %% Graphic settings
        % Apariencia de las lineas en la grafica de error
        [Perr(1).Color,Perr(2).Color,Perr(2).YColor,F1(2).Parent.YColor,F1(1).LineWidth,...
            F1(1).Color,F1(1).DisplayName,F1(2).LineWidth,F1(2).Color,F1(2).DisplayName,...
            F2.LineWidth,F2.Color,fgd1.FontSize,fgd1.Box,handles.ploterror.XLabel.String]=...
            deal('none','w','b','r',2.5,'r','Minimum misfit',1.5,'m','Current misfit',2.5,'b',10,'off','Iteration');
        % Type interpreter Latex and normalized legend
        [fgd1.Interpreter,Perr(2).YLabel.Interpreter,handles.ploterror.XLabel.Interpreter,...
            Perr(1).YLabel.Interpreter]=deal('latex');
        [Perr(1).YLabel.Units,Perr(2).YLabel.Units,handles.ploterror.XLabel.Units,...
            handles.ploterror.FontUnits]=deal('normalized');
        uistack(Perr(2),'down');
        linkaxes([Perr(1),Perr(2)],'x');
        % Load Id graphics
        LoadV('F1',F1,'F2',F2,'Perr',Perr);
    else
        [F1(1).XData,F1(2).XData,F2.XData]=deal(ER.kk);% NUmber of iterations
        [F1(1).YData,F1(2).YData,F2.YData]=deal(ER.err,ER.err1,ER.temer);% 'Minimum misfit','Current misfit','Temperature'
        linkaxes([Perr(1),Perr(2)],'x');
        [handles.ploterror.YLimMode,handles.ploterror.XLimMode]=deal('auto');
        if Tinv==2 || Tinv==3
            Perr(2).YLabel.String='Temperature';
        elseif Tinv==4 || Tinv==1 || Tinv==5
            [Perr(2).YLabel.String,Perr(2).YTickMode,Perr(2).YTick]=deal(' ','manual',[]);
        end
    end
    %Clean Memory
    clear HVT DCT mod L2i HVTp DCTp modp L2ip
    % Eliminando archivos de lectura
    delete('etc\*.txt');
    
    %% GUI graphical HV,Disp curve & models
    delete(findobj(handles.PROFILE,'type','axes'))
    ppv=subplot(1,1,1,'Parent',handles.PROFILE);
    [BTA,ALFA,~,E]=MODELO(aopt);
    G.ptlmod=plot(ppv,BTA,E,'-r',ALFA,E,':r','LineWidth',3);
    
    %% Graphic settings
    [ppv.YLim,ppv.XLabel.String,ppv.YLabel.String,ppv.YDir,ppv.FontUnits]=...
        deal([0 A.limd ],'$Velocity$ [${m} \over {s}$]','$Depth$ [$m$]','reverse','normalized');
    [ppv.YLabel.Interpreter,ppv.XLabel.Interpreter]=deal('latex');
    l3=legend(ppv,'V_s','V_p');
    l3.FontSize=8;
    %Graficando Curva HV (si existe)
    if  ~isempty(C.HVFobs)
        if length(G.Pl.HV.Children)>1
            delete(G.Pl.HV.Children(1));
        end
        G.phv=plot(G.Pl.HV,C.HVFobs,HVopt,'-r','LineWidth',4);
        l1=legend (G.Pl.HV,C.HVFile,'Best fitting HV');
        [l1.Interpreter,l1.FontSize,G.phv.Parent.YLimMode]=deal('none',8,'auto');
    end
    %Graficando Curva de dispersion (si existe)
    if  ~isempty(C.DCFobs)
        if length(G.Pl.DC.Children)>1
            delete(G.Pl.DC.Children(1));
        end
        G.pdc=plot(G.Pl.DC,C.DCFobs,DCopt,'-r','LineWidth',4);
        l1=legend (G.Pl.DC,C.DCFile,'Best fitting DC');
        [l1.Interpreter,l1.FontSize,G.pdc.Parent.YLimMode]=deal('none',8,'auto');
    end
    drawnow
end
%%  Variable initialization for iteration


%% Getting value of temperature depending on user's settings¿
% mod=zeros(A.NC,4,MM);% Modelos
% L2i2=zeros(1,MM);% Norma
% deltaE=L2i2;% Delta de energia SA
% out=zeros(A.NC,4,MM);% Modelo perturbado de salida
% DATERROR.lim=A.limd;% Limite de perfil Halfspace
% mac=0; %contador misfit
% Tfin=0;% Temperatura final
% resp=1;% Respuesta Usuario
% HVT=zeros(C.nmHV,MM);% HV salida function CFT
% DCT=zeros(C.nmDC,MM);% DC salida function  CFT
% ER.DCALL=[ER.DCALL zeros(C.nmDC,It)];% All curves DC
% ER.HVALL=[ER.HVALL zeros(C.nmHV,It)];% All curves HV
% ER.asens=cat(3,ER.asens, zeros(A.NC,4,It));% All Models
% ER.L2sens=[A.L2opt ;zeros(It,1)];% All misfit
% ER.MODALL=cat(3,ER.MODALL, zeros(A.NC,4,It));% All models depending of misfit
% ER.L2ALL=[ER.L2ALL ;zeros(It,1)];% all misfit
% ER.err=[ER.err; zeros(It,1)*nan];% Grafica de error Misfit
% ER.err1=[ER.err1; zeros(It,1)*nan];% Current misfit
% ER.temer=[ER.temer; zeros(It,1)*nan];% temperature
% ER.kk=[ER.kk ;zeros(It,1)*nan];% Iterations
% G.F1=F1;% Grafica Misfit
% G.F2=F2;% Grafica Curren misfit
[mod,L2i2,deltaE,out,DATERROR.lim,mac,Tfin,Tfin2,resp,HVT,DCT,ER.DCALL,...
    ER.HVALL,ER.asens,ER.L2sens,ER.MODALL,ER.L2ALL,ER.err,ER.err1,...
    ER.temer,ER.kk,G.F1,G.F2]=deal(zeros(A.NC,4,MM),zeros(1,MM),zeros(1,MM),zeros(A.NC,4,MM),A.limd,0,0,0,1,...
    zeros(C.nmHV,MM),zeros(C.nmDC,MM),[ER.DCALL zeros(C.nmDC,It)],[ER.HVALL zeros(C.nmHV,It)],...
    cat(3,ER.asens, zeros(A.NC,4,It)),[A.L2opt ;zeros(It,1)],cat(3,ER.MODALL, zeros(A.NC,4,It)),...
    [ER.L2ALL ;zeros(It,1)],[ER.err; nan(It,1)],[ER.err1; nan(It,1)],[ER.temer; nan(It,1)],...
    [ER.kk ;nan(It,1)],F1,F2);
% Asignando valor de temperatura (SA)
if A.T0==1
    [T1,Tini,TSAM]=deal(-A.EA*A.L2optimo/log(A.AP)); % Computing the temperature dependence from mismit
else
    [T1,Tini,TSAM]=deal(A.AP);
end

%% Kernel of SA method
handles.textmisfit.Visible= 'on';

while resp==1;
    if Tinv==1 || Tinv ==2 ||Tinv==3
        if handles.STOP.UserData==1%User stop
            break;
        end
        for j=1:A.Nsubiter+2
            if handles.STOP.UserData==1 %User stop
                break;
            end
            %% Loop for temperatures
            if j==A.Nsubiter+2
                [T1,A.aini,A.Niter,L2i1]=deal(TSAM,aopt,A.Niterf,A.L2optimo);
            end
            %% Loop for iterations
            for i = 1:A.Niter;
                if handles.STOP.UserData==1%User stop
                    break;
                end
                %% Semi-random generation of HV or Disp curves ( direct problem )
                tic;% Inicio del tiempo de cálculo para estadistica
                if A.open==0 || Tinv~=3
                    for ik=1:MM;
                        [mod(:,:,ik),HV,DC,out(:,:,ik),L2i2(:,ik)]=CTF(A,C,ik,2,G);    %%%%%%%%%%%%%
                        [DCT(:,ik),HVT(:,ik),deltaE(ik)]=deal(DC,HV,L2i2(:,ik)-L2i1);
                    end
                else
                    parfor ij=1:MM;
                        [mod(:,:,ij),HV,DC,out(:,:,ij),L2i2(:,ij)]=CTF(A,C,ij,2,G);
                        [DCT(:,ij),HVT(:,ij),deltaE(ij)]=deal(DC,HV,L2i2(:,ij)-L2i1);
                    end
                end
                Tfin=Tfin+toc; %tiempo de cálculo para estadistica
                %% Assigning variables and control parameters
                [A.flag(1,1),A.NN1,A.NN1,Cont.iter]=deal(0,A.NN2,A.NN/100,Cont.iter+1);
                [~,b2]=sort(deltaE);
                for ij=1:MM
                    % Assignin models
                    Cont.all=Cont.all+1;
                    [ER.MODALL(:,:,Cont.all),ER.L2ALL(Cont.all),ER.DCALL(:,Cont.all),...
                        ER.HVALL(:,Cont.all),ER.asens(:,:,Cont.all),ER.L2sens(Cont.all)]=...
                        deal(mod(:,:,MM),L2i2(:,MM),DCT(:,b2(1)),HVT(:,b2(1)),mod(:,:,b2(1)),L2i2(:,b2(1)));
                end
                %% Acceptance condition in SA
                modapro=ceil(rand()*length(b2));
                if deltaE(b2(1)) < 0
                    [A.aini,L2i1]=deal(mod(:,:,b2(1)),L2i2(b2(1)));
                    %% Acceptance condition in SA & MSA
                elseif rand < exp(-1*deltaE(modapro)/T1) && Tinv~=1;
                    [A.aini,L2i1]=deal(mod(:,:,b2(modapro)),L2i2(:,b2(modapro)));
                    %% Acceptance condition in Metropoli
                elseif rand < exp(-0.5*deltaE(modapro)) && Tinv==1
                    [A.aini,L2i1]=deal(mod(:,:,b2(modapro)),L2i2(:,b2(modapro)));
                end
                %% Optimal model and reallocation of control parameters
                if L2i2(:,b2(1)) < A.L2optimo
                    if Tinv==3
                        [A.flag(1,1),A.NN1]=deal(1,A.NN2/2);
                    end
                    % Asignando nuevos valores para la siguiente iteracion
                    [A.in,aopt,A.aini,A.L2optimo,L2i1,DCopt,HVopt]=deal(out(:,:,b2(1)),mod(:,:,b2(1)),mod(:,:,b2(1)),...
                        L2i2(:,b2(1)),A.L2optimo,DCT(:,b2(1)),HVT(:,b2(1)));
                    if Tinv==3 && j~=A.Nsubiter+2
                        T1=T1*A.RT;
                    end
                    %% ploting THE BEST model and curves using SA, SAM & Monte Carlo (Metropolis)
                    [BTA,ALFA,~,E]=MODELO(aopt);
                    [G.ptlmod(1).YData,G.ptlmod(2).YData,G.ptlmod(1).XData,G.ptlmod(2).XData]=...
                        deal(E,E,BTA,ALFA);
                    if  ~isempty(C.HVFobs)
                        G.phv.YData=HVopt;
                    end
                    if  ~isempty(C.DCFobs)
                        G.pdc.YData=DCopt;
                    end
                    [G.F1(1).XData,G.F1(2).XData,G.F2.XData]=deal(ER.kk);
                    [G.F1(1).YData,G.F1(2).YData,G.F2.YData]=deal(ER.err,ER.err1,ER.temer);
                    G.F1(1).Parent.XLimMode='Auto';
                    drawnow
                end
                %% Graphing of error, misfit & temperature
                if j==A.Nsubiter+2 || Tinv==1
                    mac=mac+1;% contador para la grafica del misfit
                    [DATERROR.mod(:,:,mac), DATERROR.Mis(mac),DATERROR.temp]=deal(A.aini,L2i1,TSAM);
                end
                %                 %% Graphing evolution of the error
                [ER.kk(Cont.iter),ER.err1(Cont.iter),ER.err(Cont.iter)]=deal(Cont.iter,L2i1,A.L2optimo);
                if Tinv~=1
                    ER.temer(Cont.iter)=T1;
                else
                    ER.temer(Cont.iter)=nan;
                end
                %Adjusting Control Parameters of HV curve
                if  ~isempty(C.HVFobs)
                    C.dk=C.dk+1;% Aumentando gradualmente el parametro de dk para el HV
                    if C.dk >C.mdk % dk min supera al dk maximo, se quedan igual
                        C.dk=C.mdk;
                    end
                end
                % Pintando el numero de iteraciones.
                set(handles.textmisfit,'String',['#Evaluated models','  ',num2str(Cont.all),'  ','Min. Misfit','  ',num2str(A.L2optimo)]);
                drawnow
                Tfin2=Tfin2+toc; %tiempo final de cálculo para estadistica
            end
            %% settings of temperatures
            if T1<TSAM
                TSAM=T1;
            end
            if Tinv==3 && j~=A.Nsubiter+2;
                if A.Nsubiter~=0;
                    T1=T1*1.5;
                end
            elseif Tinv==2 && j~=A.Nsubiter+2;
                if A.Nsubiter~=0;
                    T1=T1*A.RT;
                end
            end
        end
        
    else
        %% Linear Method Simplex Downhill & Interior point
        if handles.STOP.UserData==1
            break;
        end
        %% Semi-random generation of HV or Disp curves ( direct problem )
        [aopt(:,:),optimValues,BestL2f,ER,HVopt,DCopt,Cont,TfinL,Tfin2L]=Optimize(A,C,1,G,ER,handles,Tinv,Cont);
        Tfin2=Tfin2+Tfin2L; %Teimpo final de cálculo
        Tfin=TfinL;%Teimpo final de cálculo de cada modelo acomulativo
        questdlg(optimValues.message,' ','Ok','Ok');
        A.L2optimo=BestL2f; % mejor modelo encontrado
    end
    
    %% graphing of final models and curves
    [BTA,ALFA,~,E]=MODELO(aopt);
    %Graficando modelos con el mejor ajuste (Curvas HV, DC y perfiles de
    %velocidad
    [G.ptlmod(2).XData,G.ptlmod(2).YData,G.ptlmod(1).XData,G.ptlmod(1).YData]=...
        deal(ALFA,E,BTA,E);
    if  ~isempty(C.HVFobs)
        G.phv.YData=HVopt;
    end
    if  ~isempty(C.DCFobs)
        G.pdc.YData=DCopt;
    end
    handles.textmisfit.String=['# Evaluated models',' ',num2str(Cont.all),' ','Min. Misfit',' ',num2str(A.L2optimo)];
     [G.F1(1).XData,G.F1(2).XData,G.F2.XData]=deal(ER.kk);
     [G.F1(1).YData,G.F1(2).YData,G.F2.YData]=deal(ER.err,ER.err1,ER.temer);
     G.F1(1).Parent.XLimMode='Auto';
    drawnow
    
    %% Concatenation of data
    %COncatenacion de los datos listos para otro tipo de inversion.
    ind=find(ER.L2sens~=0 & ~isinf(ER.L2sens));
    Cont.iter=length(find(~isnan(ER.err) & ~isinf(ER.err)));
    [A.aini,Cont.valid,ER.HVALL,ER.DCALL,ER.asens,ER.L2sens,ER.MODALL,ER.L2ALL]=...
        deal(aopt,length(ind),ER.HVALL(:,ind),ER.DCALL(:,ind),ER.asens(:,:,ind),ER.L2sens(ind),...
        ER.MODALL(:,:,ind),ER.L2ALL(ind));
    %% Continuation of inversion with other different method
    if handles.STOP.UserData==0
        choice = questdlg('Continue Inversion?',' ','Yes','No','Yes');
        switch choice
            case 'Yes'
                resp=1;
                handles.STOP.UserData=0;
                %NN=num2str(A.NN);
                AA=selectmethod;
                uiwait(AA)
                Tinv=evalin('base','Tinv');
                % settings of type inversion (SA)
                if Tinv==2;
                    prompt={['M',char(225),'rkov`s chain length'],'Number of temperatures',['Last M',char(225),'rkov`s chain length'],'Perturbation range (%)'};
                    NUM=str2double(inputdlg(prompt,'',1,{'100','0','100','5'}));
                    if ~isempty(NUM)
                        [resp,Perr(2).YLabel.String, Perr(2).YTickMode]=deal(1,'Temperature','auto');
                        if A.T0==1
                            [Tini,T1]=deal(-A.EA*A.L2optimo/log(A.AP));
                        else
                            [T1,Tini,]=deal(A.AP);
                        end
                    else
                        [resp,A.Niter]=deal(0,100);
                    end
                elseif Tinv==0
                    [resp,A.Niter]=deal(0,100);
                    % settings of type inversion (SAM)
                elseif Tinv==3;
                    prompt={'Number of iterations','Number of Reheatings','Number of Last iterations','Perturbation range (%)'};
                    NUM=str2double(inputdlg(prompt,' ',1,{'100','0','100','5'}));
                    if ~isempty(NUM)
                        [resp,Perr(2).YLabel.String, Perr(2).YTickMode]=deal(1,'Temperature','auto');
                        if A.T0==1
                            [Tini,T1]=deal(-A.EA*A.L2optimo/log(A.AP));
                        else
                            [T1,Tini,]=deal(A.AP);
                        end
                        
                        if A.open==0
                            choice = questdlg(sprintf('The inversion method is designed to run in parallel.\n\nTo use this method in parallel?'));
                            if strcmp('Yes',choice);
                                pop=GetV('NUMW');
                                NUMW=str2double(cell2mat(inputdlg({sprintf('Enter number of Workers (Labs): (>=2 & <512)')},'',1,{num2str(pop)})));
                                if isnan(NUMW) || NUMW<2
                                else
                                    %% Parallel settings & GUI
                                    mss=msgbox({sprintf('The Parallelization Settings will take a few seconds.\n\nPlease Wait')},'HV-inv');
                                    delete(findobj(mss,'string','OK'));
                                    delete(findobj(mss,'style','frame'));
                                    set(mss, 'CloseRequestFcn','')
                                    %edit GUI
                                    [handles.Poff.Checked,handles.edit1.Enable,handles.TDC.Enable,...
                                        handles.THV.Enable,handles.FHV.Enable,handles.THVDC.Enable,...
                                        handles.PS.Enable,handles.Models.Enable,handles.Inversion.Enable,...
                                        handles.LP.Enable,handles.SPara.Enable,handles.START.Enable,...
                                        handles.STOP.Enable,handles.Menu.Enable]=deal('off');
                                    % Habilita los nucleos
                                    parallel(NUMW);
                                    % Edit GUI
                                    [handles.textNumwork.Visible,handles.Pon.Checked,handles.TDC.Enable,...
                                        handles.THV.Enable,handles.FHV.Enable,handles.THVDC.Enable,...
                                        handles.PS.Enable,handles.Models.Enable,handles.Inversion.Enable,...
                                        handles.LP.Enable,handles.SPara.Enable,handles.START.Enable,...
                                        handles.STOP.Enable,handles.Menu.Enable,handles.edit1.Enable]=deal('on');
                                    [handles.textNumwork.String,A.open,A.flag,MM]=...
                                        deal(['Connection with',' ',num2str(NUMW),' ', 'Labs.'],1,zeros(NUMW,1),NUMW);
                                    LoadV('NUMW',NUMW,'open',1);
                                    delete(mss);
                                end
                            end
                        else
                            MM=GetV('NUMW');
                            A.flag=zeros(MM,1);
                        end
                    else
                        [resp,A.Niter]=deal(0,100);
                    end
                    % settings of type inversion (Metropoli)
                elseif Tinv==1;
                    NUM=str2double(inputdlg({'Number of iterations','Perturbation range (%)'},'',1,{'100','5'}));
                    if ~isempty(NUM)
                        NUM(4)=NUM(2);
                        [resp,NUM(2),NUM(3),Perr(2).YLabel.String,Perr(2).YTickMode,Perr(2).YTick]=...
                            deal(1,0,0,' ','manual',[]);
                    else
                        [resp,A.Niter]=deal(0,100);
                    end
                    % settings of type inversion (Simplex Downhill)
                elseif Tinv==4;
                    prompt={'Maximum number of iterations allowed','Maximum number of evaluations allowed','Termination tolerance on the function value (TolFun)','Termination tolerance on parameters (TolX)'};
                    NUM=str2double(inputdlg(prompt,' ',1,{'100','200','1e-4','1e-4'}));
                    if ~isempty(NUM)
                        [resp,Perr(2).YLabel.String,Perr(2).YTickMode,Perr(2).YTick]=deal(1,' ','manual',[]);
                    else
                        [resp,A.Niter]=deal(0,100);
                    end
                    % settings of type inversion (Interior point)
                elseif Tinv==5;
                    prompt={'Maximum number of iterations allowed','Maximum number of evaluations allowed','Termination tolerance on the function value (TolFun)','Termination tolerance on parameters (TolX)'};
                    NUM=str2double(inputdlg(prompt,',',1,{'100','200','1e-4','1e-4'}));
                    if ~isempty(NUM)
                        [resp,Perr(2).YLabel.String,Perr(2).YTickMode,Perr(2).YTick]=deal(1,' ','manual',[]);
                    else
                        [resp,A.Niter]=deal(0,100);
                    end
                end
            otherwise
                [resp,A.Niter]=deal(0,100);
        end
        %% reassignment of settings for inversion type
        if resp==1
            [A.aini,A.Niter, A.Nsubiter,A.Niterf,A.NN]=deal(aopt,NUM(1),NUM(2),NUM(3),NUM(4));
            %% concatenating the curve data & information from error
            if Tinv==1 || Tinv==2 || Tinv==3
                [It,ER.kk,ER.err,ER.err1,ER.temer,ER.DCALL,ER.HVALL,ER.asens,ER.L2sens,...
                    ER.MODALL,ER.L2ALL,A.ST]=deal(It+((A.Nsubiter+1)*A.Niter+A.Niterf),...
                    [ER.kk; nan((A.Nsubiter+1)*A.Niter+A.Niterf,1)],[ER.err; nan((A.Nsubiter+1)*A.Niter+A.Niterf,1)],[ER.err1; nan((A.Nsubiter+1)*A.Niter+A.Niterf,1)],...
                    [ER.temer; nan((A.Nsubiter+1)*A.Niter+A.Niterf,1)],[ER.DCALL zeros(C.nmDC,(A.Nsubiter+1)*A.Niter+A.Niterf)],...
                    [ER.HVALL zeros(C.nmHV,(A.Nsubiter+1)*A.Niter+A.Niterf)],cat(3,ER.asens, zeros(A.NC,4,(A.Nsubiter+1)*A.Niter+A.Niterf)),...
                    [ER.L2sens ;zeros((A.Nsubiter+1)*A.Niter+A.Niterf,1)],cat(3,ER.MODALL, zeros(A.NC,4,(A.Nsubiter+1)*A.Niter+A.Niterf)),...
                    [ER.L2ALL ;zeros((A.Nsubiter+1)*A.Niter+A.Niterf,1)],A.ST+A.Nsubiter);
            end
            %Remove Statistical models
            if isfield(DATERROR, 'mod')
                DATERROR = rmfield(DATERROR,'mod');
                DATERROR = rmfield(DATERROR,'Mis');
                DATERROR = rmfield(DATERROR,'temp');
            end
            mac=0;% Inicia contador del misfit
        else
            %HAbilitando botones de la GUI
            [handles.STOP.Visible,handles.smodel.Visible,handles.START.Visible]= deal('off','on','on');
            % Mensage informativo del misfit minimo
            msgbox(['Finished program. ',' The minimun misfit is ' num2str(A.L2optimo)],' ');
        end
    end
end
%% Termiantion of inversion
% Assigning variables in the workspace
if ~isempty(ini1)
    [ER.MODALL(:,:,end+1),ER.asens(:,:,end+1),ER.HVALL(:,end+1),ER.DCALL(:,end+1),...
        ER.L2ALL(end+1),ER.L2sens(end+1) ]=deal(modini1,modini1,HVTini1,DCTini1,L2ini1,L2ini1);
end
A.aini=aopt;
% Assigning variables in the workspace
LoadV('aopt',aopt,'A',A,'HVALL',ER.HVALL,'DCALL',ER.DCALL,'asens',ER.asens,'L2sens',ER.L2sens,...
    'MODALL',ER.MODALL,'L2ALL',ER.L2ALL,'DATAERROR',DATERROR);
%% Report from final Inversion
if Tinv==3
    ST1='Modified SA';
    KGF=['Number of reheatings : ',num2str(A.ST)];
elseif Tinv==2
    ST1='SA';
    KGF=['Number of temperatures : ',num2str(A.ST)];
elseif Tinv==1
    %Reporte final para Metropoli Walk
    if pass==1
        ST1='Metropoli Walk';
        handles.Report.Enable= 'on';
        handles.Report.String=[{'Inversion report'};{' '};{['Method Inversion : ',ST1 ]};{['Number of iterations : ',num2str(Cont.iter)]};{'  '};...
            {['Perturbation range : ',num2str(A.NN),'%']};{' '};...
            {['Total models : ',num2str(Cont.all)]};{['Accepted models  : ',num2str(Cont.valid)]} ;{['Models for statistics : ',num2str(mac)]};{'  '};...
            {['Maximum misfit : ',num2str(max(ER.L2sens))]};{['Minimum misfit : ',num2str(min(ER.L2sens))]};{'  '};...
            {['Total time : ',num2str(Tfin2)]};{['Avg time per model : ',num2str(Tfin/Cont.all)]}];
    end
end
if Tinv==2   || Tinv==3
    if pass==1
        %Reporte final para SA y SAM
        handles.Report.Enable= 'on';
        handles.Report.String=[{'Inversion report'};{' '};{['Method Inversion : ',ST1 ]};{['Number of Iterations : ' ,num2str(Cont.iter)]};{KGF};{['Number of last iterations : ',num2str(mac)]};{'  '};...
            {['Initial Temperature : ',num2str(Tini)]};{['Final Temperature : ',num2str(T1)]};{['Temperature ratio : ', num2str(A.RT)]} ;{['Perturbation range : ',num2str(A.NN),'%']};{' '};...
            {['Models evaluated : ',num2str(Cont.all)]};{['Valid models  : ',num2str(Cont.valid)]} ;{['Models for statistics : ',num2str(mac)]};{'  '};...
            {['Maximum misfit : ',num2str(max(ER.L2sens))]};{['Minimum misfit : ',num2str(min(ER.L2sens))]};{'  '};...
            {['Total time : ',num2str(Tfin2)]};{['Avg time per model : ',num2str(Tfin/Cont.all)]}];
    end
end
%Reporte final para los metodos de inversion lineal.
if Tinv==4
    ST1='Downhill simplex ';
elseif Tinv==5
    ST1='Interior-point';
end
if pass==1 && (Tinv==4 || Tinv==5);
    handles.Report.Enable= 'on';
    handles.Report.String=[{'Inversion report'};{' '};{['Method Inversion : ',ST1 ]};{['Number of iterations : ',num2str(Cont.iter)]};{'  '};...
        {['Valid models : ',num2str(Cont.valid)]};{['Models evaluated  : ',num2str(Cont.all)]} ;{'  '};...
        {['Maximum misfit : ',num2str(max(ER.L2sens))]};{['Minimum misfit : ',num2str(min(ER.L2sens))]};{'  '};...
        {['Total time : ',num2str(Tfin2)]};{['Avg time per model : ',num2str(Tfin/Cont.all)]}];
end
