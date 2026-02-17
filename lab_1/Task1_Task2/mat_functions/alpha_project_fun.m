function [alpha_project] = alpha_project_fun(x_camber,y_camber,chord)
% This function computes the theorical project angle of the profile

% First change of variables of the integral  [0;l] --> [-l/2;+l/2]
s_camber = x_camber - (chord/2);

% Derivate y'
dy_ds = gradient(y_camber, s_camber);

% Function
g = @(s) interp1(s_camber,dy_ds,s,'pchip');

% Grid of eta
N = 2000;
eta = linspace(0,pi,N);

% Second change of variables [s = -l/2 * cos(eta)]
s_eta = -(chord/2).*cos(eta);

% Computing the integral:  alpha_project = 1/pi * Int[y'(s(eta))][0;pi]
alpha_project = (1/pi) * trapz(eta,g(s_eta));


end