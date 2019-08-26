function [P] = fitthreemodelprocess(data,t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define Two model process function 
fh = @(x,p) p(1)*(1-exp(-x/p(2)))+p(3)*(1-exp(-x/p(4)))+p(5)*(1-exp(-x/p(6)));
% p = [ Af, Tf, As, Ts, Am, Tm]

%Define Error Function
errfh = @(p,x,y) sum((y(:)-fh(x(:),p)).^2);

% Initial Guess 
p0 = [ 0.9, 3, 2, 20, 1.2, 15];

% search for solution - minimization of error function 
P = fminsearch(errfh,p0,[],t,data);

% plot results 
% figure
% plot(t,data,'b.',t,fh(t,P),'r-');

end
