%% Vs_30
% compute the Vs30 of a velocity profile

    % INPUTS
    % d - a vector of depths
    
    % v - a  vector of velocities (same length as d)
    
    % OUTPUTS
    
    % vs_30 - Vs_30 of the profile
    
    % class - NEHRP site classification of the profile
    
    
%% Author: Marshall Pontrelli
% Date: 12/7/2020

function [vs_30, class] =  Vs_30(d,v)
    dd = 0; % depth counter to keep track of depth
    counter = 0;
    if sum(d) >= 30
        while dd < 30
            counter = counter + 1;
            dd = dd + d(counter);
        end
        for i = 1:counter - 1
            den(i) = d(i)/v(i);
        end
        den(counter) = (30 - sum(d(1:counter-1)))/v(counter);
        vs_30 = 30/sum(den);
       a = strcat('Vs30 =',{' '},num2str(vs_30),{' '},'m/s');
       a = a{1};
       %disp(a)
       if vs_30 > 1500
           class = 'A';
       end
       if vs_30 > 760 &&  vs_30 <=1500
           class = 'B';
       end
       if vs_30 > 360 &&  vs_30 <=760
           class = 'C';
       end
       if vs_30 > 180 &&  vs_30 <=360
           class = 'D';
       end
       if vs_30 <= 180
           class = 'E';
       end
       b = strcat('class =',{' '},class);
       b = b{1};
       %disp(b)
       
       
    else
        vs_30 = [];
        disp('Profile is less than 30 meters')%
    end
end