function GEidx = get_GEidx(time, GE)
% Finds the indices for the gait events for a given time vector and gait
% event times, in GE.
%
% INPUT
% GE, gait event times. Assumes each row is 1 set of gait events and the
% columns are the different gait events.
%
% time, time vector for which to identify the gait event indices
% Assumes time(1) is the same actual time. i.e. synced.
%
% OUTPUT
% GEidx is a matrix of the same size as GE.
% GEidx is the indices that match the times in GE.

% preallocate;
GEidx = nan(size(GE));

if ~isempty(find(time == GE(1,1), 1))
    
    % can find exact times. Ex. time and GE from the same system
    for ev = 1:length(GEidx)
        GEidx(:,ev) = find(ismember(time, GE(:,ev)));
    end
else
    % can't find exact times. Ex. time for EMG and GE from Mocap
    for ev = 1:size(GE,2)
        for stridenum = 1:size(GE,1)
            GEidx(stridenum, ev) = find(time >= GE(stridenum,ev),1,'first');
        end
    end   
end
