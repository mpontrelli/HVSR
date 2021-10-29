
% Compute and show mean model, std deviation, covariance matrix and 
% and correlation matrix
%
%     Copyright (C) 2014,2016 Antonio García-Jerez, José Piña-Flores.
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
%

DATAERROR=GetV('DATAERROR'); % Load DATAERROR structure with a list of models and misfits
if isfield(DATAERROR, 'Mis')
    
    % Settings for representation
    linewidth=5;
    meanP='k';col_errbar_P=[0 0 0];
    meanS='k';col_errbar_S=[0 0 0];

    % Take number of models and layers
    ndatos=length(DATAERROR.Mis); % number of models
    ncapas=size(DATAERROR.mod(:,:,1),1);% number of layers, including halfspace

    % Select model parameters for statistics (from left side to right side and from up to down)
    icoords=[1:4:4*(ncapas-1) 2:4:4*ncapas 3:4:4*ncapas ]; % all thicknesses and P and S wave velocities
    simbolos{1}='h';simbolos{2}='Vp';simbolos{3}='Vs';simbolos{4}='Dens';% for labels
    dim=length(icoords);

    % Prepare matrix "coords" for statistics (each row is a model, each column is a selected model parameter)
    coords=nan(ndatos,length(icoords));
    for index=1:ndatos
        modelo=DATAERROR.mod(:,:,index).';
        coords(index,:)=modelo(icoords);
    end
    
    % Find best model and compute mean model
    [kk ,pos]=min(DATAERROR.Mis);
    best_model=DATAERROR.mod(:,:,pos).';
    modelo_promedio=mean(coords,1);

    % Covariance matrix (matriz_cov)
    matriz_cov=zeros(dim);
    for index=1:ndatos,
        matriz_cov=matriz_cov+transpose(coords(index,:)-modelo_promedio)*(coords(index,:)-modelo_promedio);
    end
    matriz_cov=matriz_cov/ndatos;
    
    % Show the model space if it is two-dimensional
%     if dim==2,
%         figure;axis;hold on;box on;
%         plot(coords(:,1),coords(:,2),'k.')
%         plot(modelo_promedio(1),modelo_promedio(2),'bo')
%         plot(best_model(icoords(1)),best_model(icoords(2)),'r*')
%         xlabel([simbolos{rem(icoords(1)-1,4)+1},num2str(ceil(icoords(1)/4))])
%         ylabel([simbolos{rem(icoords(2)-1,4)+1},num2str(ceil(icoords(2)/4))])
%         oid=line(modelo_promedio(1)*[1 1],modelo_promedio(2)+sqrt(matriz_cov(2,2))*[-1 1]);
%         set(oid,'Color',[0 0 0]);
%         oid=line(modelo_promedio(1)+sqrt(matriz_cov(1,1))*[-1 1],modelo_promedio(2)*[1 1]);
%         set(oid,'Color',[0 0 0]);    
%         axis square
%     end

    % Compute and show correlation matrix
    matriz_corrs=matriz_cov(1:dim+1:end);% take diagonal from matriz_cov as a row [C11,C22,...]
    matriz_corrs=sqrt(matriz_corrs.' * matriz_corrs);% transform into matrix sqrt(Cii*Cjj)
    matriz_corrs=abs(matriz_cov)./matriz_corrs;% correlation matrix
    matriz_corrs_ampliada=zeros(size(matriz_corrs)+1);% a trick that avoid missing a column and a row in the representation
    matriz_corrs_ampliada(1:end-1,1:end-1) = matriz_corrs;
    matriz_corrs_ampliada(end,1:end-1)     = matriz_corrs_ampliada(end-1,1:end-1);
    matriz_corrs_ampliada(1:end-1,end)     = matriz_corrs_ampliada(1:end-1,end-1);
    matriz_corrs_ampliada(end,end)         = matriz_corrs_ampliada(end,end-1);
    figure;eje=axes;surf(matriz_corrs_ampliada);view([0 -90]);
    tickcell={dim};
    for index=1:dim,tickcell{index}=[simbolos{rem(icoords(index)-1,4)+1},num2str(ceil(icoords(index)/4))];end
    axis equal;axis tight;caxis([0 1]);colormap(1-gray);
    colorbar;
    set(eje,'xtickmode','manual','xticklabelmode','manual','xtick',(1:dim)+0.5,'xticklabel',tickcell);
    set(eje,'ytickmode','manual','yticklabelmode','manual','ytick',(1:dim)+0.5,'yticklabel',tickcell,'fontsize',14);
    set(eje,'ytickmode','manual','yticklabelmode','manual','ytick',(1:dim)+0.5,'yticklabel',tickcell,'fontsize',14);
    set(eje,'XAxisLocation','top');
    
    % Plot mean Vs profile and std deviation
    modelo_medio_completo(icoords)=modelo_promedio;
    d    = modelo_medio_completo(1:4:(ncapas-1)*4);
    bta  = modelo_medio_completo(3:4:ncapas*4);
    figure(5);subplot(1,3,2);hold on
    hg=evalin('base','hg');
    DATAERROR=evalin('base','DATAERROR');
    y(2:2:2*length(d)+2)=[cumsum(d) DATAERROR.lim];
    y(3:2:2*length(d)+1)=y(2:2:2*length(d));
    x(2:2:2*length(d)+2)=bta;
    x(1:2:2*length(d)+1)=bta;
    dibujo=plot(x,y,'color',[.0 .0 .0],'LineWidth',linewidth);
    % Axes, limits and legend
    set(gca,'YDir','reverse');
    ylim([0 y(end)]);
    l4=legend([hg(2) dibujo],'Best model V_s','Mean model V_s','Location','best');
    set(l4,'FontSize',10);    
    for index=1:length(icoords),
        kk=find(icoords(index)==1:4:(ncapas-1)*4);
        if ~isempty(kk)
            oid=line([1 1]*mean(bta(kk:kk+1)),[-1 1]*sqrt(matriz_cov(index,index))+sum(d(1:kk))); % plot standard deviation of layer thickness
            set(oid,'linewidth',2,'color',col_errbar_S);% set color col_errbar_S
            continue;        
        end
        kk=find(icoords(index)==(3:4:ncapas*4));
        if ~isempty(kk)
            if kk<ncapas
                oid=line([-1 1]*sqrt(matriz_cov(index,index))+bta(kk),[1 1]*(sum(d(1:kk-1))+d(kk)/2));% plot standard deviation of Vs
            else
                oid=line([-1 1]*sqrt(matriz_cov(index,index))+bta(kk),[1 1]*(sum(d(1:kk-1))+(max(y)-sum(d(1:kk-1)))/2));% Same for halfspace            
            end
            set(oid,'linewidth',2,'color',col_errbar_S);% set color col_errbar_S
            continue;        
        end
    end

    % Plot mean Vp profile and std deviation
    alfa = modelo_medio_completo(2:4:ncapas*4);
    figure(5);subplot(1,3,1);hold on;
    x(2:2:2*length(d)+2)=alfa;
    x(1:2:2*length(d)+1)=alfa;
    dibujo=plot(x,y,'color',[.0 .0 .0],'LineWidth',linewidth);    
    % Axes, limits and legend
    set(gca,'YDir','reverse');
    ylim([0 y(end)]);
    l4=legend([hg(1) dibujo],'Best model V_p','Mean model V_p','Location','best');
    set(l4,'FontSize',10)
    for index=1:length(icoords),
        kk=find(icoords(index)==1:4:(ncapas-1)*4);
        if ~isempty(kk)
            oid=line([1 1]*mean(alfa(kk:kk+1)),[-1 1]*sqrt(matriz_cov(index,index))+sum(d(1:kk))); % plot standard deviation of layer thickness
            set(oid,'linewidth',2,'color',col_errbar_P);% set color col_errbar_P
            continue;        
        end
        kk=find(icoords(index)==(2:4:ncapas*4));
        if ~isempty(kk)
            if kk<ncapas
                oid=line([-1 1]*sqrt(matriz_cov(index,index))+alfa(kk),[1 1]*(sum(d(1:kk-1))+d(kk)/2)); % plot standard deviation of layer Vs
            else
                oid=line([-1 1]*sqrt(matriz_cov(index,index))+alfa(kk),[1 1]*(sum(d(1:kk-1))+(max(y)-sum(d(1:kk-1)))/2));% Same for halfspace            
            end
            set(oid,'linewidth',2,'color',col_errbar_P);% set color col_errbar_P
            continue;        
        end
    end

    figure(5);
    subplot(1,3,3);
    ylim([0 DATAERROR.lim]);%set depth limit for density subplot
end


