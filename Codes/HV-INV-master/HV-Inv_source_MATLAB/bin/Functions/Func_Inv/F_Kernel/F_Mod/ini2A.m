%     Copyright (C) 2014,2016 Antonio Garc�a-Jerez, Jos� Pi�a-Flores.
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

function [aale,out]=ini2A(aini,NC,var,NN,ZLVs,ZLVp,flag,in) 
%% Random generation of a model (aale) around a given one (aini) considering a perturbation range NN (maximum)
%% Ground parameter ranges (var) and constraints (ZLVs,ZLVp) are fulfilled

% Ckeck validity of the initial model
if find(aini(:,3)<var(:,7)),fprintf(1,'INI2: Vs in aini is below the minimum');end
if find(aini(:,3)>var(:,8)),fprintf(1,'INI2: Vs in aini exceeds the maximum');end
if find(aini(:,2)<var(:,4)),fprintf(1,'INI2: Vp in aini is below the minimum');end
if find(aini(:,2)>var(:,5)),fprintf(1,'INI2: Vp in aini exceeds the maximum');end    
if ~isempty(find(diff(var(:,7))<0,1))&&~ZLVs,errordlg('INI2: The lower bounds of Vs do not increase downwards in var');end
if ~isempty(find(diff(var(:,8))<0,1))&&~ZLVs,errordlg('INI2: The upper bounds of Vs do not increase downwards in var');end  
if ~isempty(find(diff(aini(:,3))<0,1))&&~ZLVs,errordlg('INI2: Vs in aini does not increase downwards');end
if ~isempty(find(diff(var(:,4))<0,1))&&~ZLVp,errordlg('INI2: The lower bounds of Vp do not increase downwards in var');end
if ~isempty(find(diff(var(:,5))<0,1))&&~ZLVp,errordlg('INI2: The upper bounds of Vp do not increase downwards in var');end  
if ~isempty(find(diff(aini(:,2))<0,1))&&~ZLVp,errordlg('INI2: Vp does not increase downwards in aini');end
%if find(aini(:,3)./aini(:,2)<sqrt((1-(2*var(:,14)))./(2*(1-var(:,14))))&aini(:,2)./aini(:,3)>sqrt((2*(1-var(:,14)))./(1-(2*var(:,14))))),
epsil=1e-3;
if find(sqrt((1-(2*var(:,14)))./(2*(1-var(:,14))))-aini(:,3)./aini(:,2)>epsil),
    fprintf(1,'INI2 :  maximum Poisson ratio exceeded in aini');
end
%if find(aini(:,3)./aini(:,2)>sqrt((1-(2*var(:,13)))./(2*(1-var(:,13))))&aini(:,2)./aini(:,3)<sqrt((2*(1-var(:,13)))./(1-(2*var(:,13))))),
if find(aini(:,3)./aini(:,2)-sqrt((1-(2*var(:,13)))./(2*(1-var(:,13))))>epsil),
    %aini(:,3)./aini(:,2)-sqrt((1-(2*var(:,13)))./(2*(1-var(:,13))))
    %sqrt((2*(1-var(:,13)))./(1-(2*var(:,13))))-aini(:,2)./aini(:,3)
    fprintf(1,'INI2 : Minimum Poisson ratio is not reached in aini');
end

% The Poisson's ratio is not limited by the perturvation range NN. Only user ranges (table) are considered. 

if ~flag
    limitedvar=[1 2 3 4];% thickness VP VS density
    var(:,limitedvar*3-2)=max(var(:,limitedvar*3-2),aini(:,limitedvar)-var(:,limitedvar*3)*NN); % Minima: fulfill both user's constraints (var) and perturbation range
    var(:,limitedvar*3-1)=min(var(:,limitedvar*3-1),aini(:,limitedvar)+var(:,limitedvar*3)*NN); % Maxima: fulfill both user's constraints (var) and perturbation range 
    var(:,limitedvar*3)=var(:,limitedvar*3-1)-var(:,limitedvar*3-2);% Update amplitudes of ranges
else % MSA
    limitedvar=[1 2 3 4];% thickness VP VS density
    for index=1:length(limitedvar)
        for kl=1:NC
            if in(kl,limitedvar(index))>=0 % The parameter have to be increased
                var(kl,limitedvar(index)*3-2)=aini(kl,limitedvar(index));% The lower bound is limited by the current model
                var(kl,limitedvar(index)*3-1)=min(var(kl,limitedvar(index)*3-1),aini(kl,limitedvar(index))+var(kl,limitedvar(index)*3)*NN); % The upper bound comes from (aini,NN) or from the global limit (var table)
            else % The parameter have to be decreased
                var(kl,limitedvar(index)*3-2)=max(var(kl,limitedvar(index)*3-2),aini(kl,limitedvar(index))-var(kl,limitedvar(index)*3)*NN); % The lower bound comes from (aini,NN) or from the global limit (var table)
                var(kl,limitedvar(index)*3-1)=aini(kl,limitedvar(index)); % The upper bound is limited by the current model                
            end                
        end
    end
    weaklimitedvar=setdiff([1 2 3 4],limitedvar);
    var(:,weaklimitedvar*3-2)=max(var(:,weaklimitedvar*3-2),aini(:,weaklimitedvar)-var(:,weaklimitedvar*3)*NN); % Minima: fulfill both user's constraints (var) and perturbation range
    var(:,weaklimitedvar*3-1)=min(var(:,weaklimitedvar*3-1),aini(:,weaklimitedvar)+var(:,weaklimitedvar*3)*NN); % Maxima: fulfill both user's constraints (var) and perturbation range 
    var(:,weaklimitedvar*3)=var(:,weaklimitedvar*3-1)-var(:,weaklimitedvar*3-2);% Update amplitudes of ranges   
    var(:,3:3:15)=var(:,2:3:14)-var(:,1:3:13);% Update amplitudes of ranges
end

all_correct=0;
not_colpased=true(1,NC);
while ~all_correct

    if ~ZLVs % Work with Vs
        % Replace var(kl,8) with the minimum upper Vs bound taken among the kl-th and deeper layers
        for kl=NC-1:-1:1
            var(kl,8)=min(var(kl,8),var(kl+1,8));
        end
        % Replace var(kl,7) with the maximum lower Vs bound taken among the kl-th and shallower layers
        for kl=2:NC
            var(kl,7)=max(var(kl,7),var(kl-1,7));
        end
    elseif ~ZLVp % Work with Vp
        for kl=NC-1:-1:1
            var(kl,5)=min(var(kl,5),var(kl+1,5));
        end
        for kl=2:NC
            var(kl,4)=max(var(kl,4),var(kl-1,4));
        end        
    end
    
    % Check again that there are solutions for the independent variable fulfilling the Poisson's ratios
    if ZLVp % Work in Vs
        
        var(not_colpased,7)=max(var(not_colpased,7),var(not_colpased,4).*sqrt((1-(2*var(not_colpased,14)))./(2*(1-var(not_colpased,14))))); % Increments minimum Vs if necessary (from Vp and Poisson ranges)
        var(not_colpased,8)=min(var(not_colpased,8),var(not_colpased,5).*sqrt((1-(2*var(not_colpased,13)))./(2*(1-var(not_colpased,13))))); % Decrements maximum Vs if necessary (from Vp and Poisson ranges)
        var(not_colpased,9)=var(not_colpased,8)-var(not_colpased,7);
        
        % Fix Vs if the amplitude of the Vs interval colapsed 
        epsil=1e0;
        not_colpased=~(var(:,9)<epsil);var(~not_colpased,7)=aini(~not_colpased,3);var(~not_colpased,8)=aini(~not_colpased,3);var(~not_colpased,9)=0;
        
        % Loop if (minVS - maxVS)>tiny
        epsil=1e-3; 
        if find(max(var(:,7),var(:,4).*sqrt((1-(2*var(:,14)))./(2*(1-var(:,14)))))-min(var(:,8),var(:,5).*sqrt((1-(2*var(:,13)))./(2*(1-var(:,13)))))>epsil)
            continue
        end        
    else % Work in Vp
        
        var(not_colpased,4)=max(var(not_colpased,4),var(not_colpased,7).*sqrt((2*(1-var(not_colpased,13)))./(1-(2*var(not_colpased,13)))));
        var(not_colpased,5)=min(var(not_colpased,5),var(not_colpased,8).*sqrt((2*(1-var(not_colpased,14)))./(1-(2*var(not_colpased,14)))));
        var(not_colpased,6)=var(not_colpased,5)-var(not_colpased,4);
        
        % Fix Vp if the amplitude of the Vp interval colapsed
        epsil=1e0;
        not_colpased=~(var(:,6)<epsil);var(~not_colpased,4)=aini(~not_colpased,2);var(~not_colpased,5)=aini(~not_colpased,2);var(~not_colpased,6)=0;
        
        % Loop if (minVP - maxVP)>tiny
        epsil=1e-3;
        if find(max(var(:,4),var(:,7).*sqrt((2*(1-var(:,13)))./(1-(2*var(:,13)))))-min(var(:,5),var(:,8).*sqrt((2*(1-var(:,14)))./(1-(2*var(:,14)))))>epsil)
            continue
        end                
    end
     
    if ~ZLVs
        % Check increasing Vs bounds        
        if ~isempty(find(diff(var(:,7))<0,1)),
            continue;
        elseif ~isempty(find(diff(var(:,8))<0,1))
            continue;
        end
    end    
    
    if ~ZLVp
        % Check increasing Vp bounds
        if ~isempty(find(diff(var(:,4))<0,1)),
            continue;
        elseif ~isempty(find(diff(var(:,5))<0,1))
            continue;
        end
    end    

    all_correct=1; % all tests passed    
end

aale=ini1_plus(NC,var,ZLVs,ZLVp);% generate a random model in the new ranges

%% Record the variation directions of the parameters. Used for MSA (Pi�a's experimental Modified version of Simulated Annealing)
if flag==0
    out=sign(aale-aini);% compute
else
    out=in; % Those in the input "in" have been preserved
end

return