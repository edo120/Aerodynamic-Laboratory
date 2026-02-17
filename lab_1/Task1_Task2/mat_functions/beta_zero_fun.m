function [beta_zero] = beta_zero_fun(x_camber,y_camber,Chord)
% This function computes beta zero angle of the profile

% First change of variables of the integral  [0;l] --> [-l/2;+l/2]
s_camber = x_camber - (Chord/2);

% Derivate y'
dy_ds = gradient(y_camber, s_camber);

% Function
g = @(s) interp1(s_camber,dy_ds,s,'pchip');

% Grid of eta
N = 2000;
eta = linspace(0,pi,N);

% Second change of variables  [s = -l/2 * cos(eta)]
s_eta = -(Chord/2).*cos(eta);

% Computing the integral:  
% beta_zero = 2/pi * Int[(sin(eta)^2)*y'(s(eta))][0;pi]
beta_zero = (2/pi) * trapz(eta,(sin(eta).^2).*g(s_eta));

end