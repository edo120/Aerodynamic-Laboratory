function [surface,camber] = camber_fun(geo,Chord,NPanels)
% camber_fun separates upper and lower surfaces curves and
% computes the camber line

% Separating upper surface curve
[~,LE] = min(geo.x);
surface.upper.y_upper = geo.y(LE:end);
surface.upper.x_upper = geo.x(LE:end);

[x_u,iu] = unique(surface.upper.x_upper);
y_u = surface.upper.y_upper(iu);

% separating lower surface curve
surface.lower.y_lower = geo.y(1:LE);
surface.lower.x_lower = geo.x(1:LE);

[x_l,il] = unique(surface.lower.x_lower);
y_l = surface.lower.y_lower(il);

% Interpolating on the same axis the curves
x_camber = linspace(0,Chord,NPanels)';
y_upper_interp = interp1(x_u,y_u,x_camber,'pchip');
y_lower_interp = interp1(x_l,y_l,x_camber,'pchip');

y_camber = (y_upper_interp + y_lower_interp)/2;
y_camber(1) = 0;
y_camber(end) = 0;

camber.x = x_camber;
camber.y = y_camber;

end