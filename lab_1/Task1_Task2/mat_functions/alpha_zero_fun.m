function [alpha_zero] = alpha_zero_fun(x_camber,y_camber,chord,alpha_project)
% This function computes the zero lift angle of the profile

% First change of variables of the integral  [0;l] --> [-l/2;+l/2]
s_camber = x_camber - (chord/2);

% Derivate y'
dy_ds = gradient(y_camber, s_camber);

% Function
g = @(s) interp1(s_camber,dy_ds,s,'pchip');

% Grid of eta
N = 2000;
eta = linspace(0,pi,N);

% Second change of variables  [s = -l/2 * cos(eta)]
s_eta = -(chord/2).*cos(eta);

% Computing the integral:  
% alpha_zero = alpha_project -1/pi * Int[cos(eta)*y'(s(eta))][0;pi]
alpha_zero = alpha_project - ( (1/pi) * trapz(eta,cos(eta).*g(s_eta)));

end