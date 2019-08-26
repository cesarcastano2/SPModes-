function [cv] = movcv(y)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% cov = movmean(y,[0 99]);

cv = movstd(y,[0 14])./movmean(y,[0 14]);

end

