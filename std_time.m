function [std_steadystate_everysix,std_steadystate_avg] = std_time(y)
size_loop = length(y);
N = 6;   % number for the size of the standard deviation it calculates 
% std_steadystate = zeros(1,size_loop);
for i=1:1:(size_loop-6)
    std_steadystate_everysix(i,1) = std(y(i:i+N));
end

std_steadystate_avg = std(y(:));