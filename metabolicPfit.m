function [P,fh] = metabolicPfit(mp,v)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Define function to fit
fh = @(x,p) p(1)*(x.^p(2))+p(3);
% p = [ Af, Tf, As, Ts]

%Define Error Function
errfh = @(p,x,y) sum((y(:)-fh(x(:),p)).^2);

% Initial Guess 
p0 = [ 0.5, 2, 1];

% search for solution - minimization of error function 
P = fminsearch(errfh,p0,[],v,mp);

% plot results 
% figure
% plot(v,mp,'b.',v,fh(v,P),'r-');

end
