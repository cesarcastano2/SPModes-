function [P,fh] = fitsinglemodelprocess_sw(data,t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define Two model process function 
% fh = @(x,p) (p(1)*exp(-(x)/p(2)))+p(3);
fh = @(x,p) (p(1)* x) + p(2);
% p = [ Af, Tf, As, Ts]

%Define Error Function
errfh = @(p,x,y) sum((y(:)-fh(x(:),p)).^2);

% Initial Guess 
p0 = [ 0.5, 0.75];

% search for solution - minimization of error function 
P = fminsearch(errfh,p0,[],t,data);

% plot results 
% figure
% plot(t,data,'b.',t,fh(t,P),'r-');

end
