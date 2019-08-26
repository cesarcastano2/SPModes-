function [Opp] = coordinate_transform(Hyp,degrees)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for i=1:length(Hyp)
Opp(i) = sin(degrees)*Hyp(i);
end
end

