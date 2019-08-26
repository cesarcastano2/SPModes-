function [P,lm] = linearreg(data,t)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define Two model process function 
fh = @(x,p) p(1)*(1-exp(-x/p(2)));
% p = [ Af, Tf, As, Ts]

%Define Error Function
errfh = @(p,x,y) sum((y(:)-fh(x(:),p)).^2);

% Initial Guess 
p0 = [ 0.9, 0.2];

% search for solution - minimization of error function 
P = fminsearch(errfh,p0,[],t,data);

tbl = table(t,P(1),P(2),data,'VariableNames',{'t','P','P','data'});
lm = fitlm(tbl,'data~p(1)*(1-exp(-t/p(2)))');

% plot results 
% figure
% plot(t,data,'b.',t,fh(t,P),'r-');

end