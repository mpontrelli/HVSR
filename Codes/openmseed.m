%% open mseed file
% rdmseed will create a structure of blockettes. This code turns these
% blockettes into vectors NS, EW and V.

    % INPUTS
    
    % x is the data structure from the rdmseed code

function [NS, EW, V] = openmseed(x)
    V = [];
    EW = [];
    NS = [];
    for i = 1:length(x)
        if x(i).ChannelIdentifier == 'ENZ'
            V = vertcat(V, x(i).d);
        end
        if x(i).ChannelIdentifier == 'ENN'
            NS = vertcat(NS, x(i).d);
        end
        if x(i).ChannelIdentifier == 'ENE'
            EW = vertcat(EW, x(i).d);
        end
    end

end