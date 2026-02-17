function [Cm] = Cm_fun(geo,CP_vect,NPanels)
% This function calculates the Cm (lower case l --> 2D airfoil) of an
% airfoil from the Leading Edge (LE) 

addpath mat_functions

% Compute elements from the panels 
[centers,normals,~,~,~,~,lengths,~,~] = CreatePanels(geo);

% Compute Cm 
Cm = 0;
c = max(geo.x) - min(geo.x);             % Chord length 

for i = 1:NPanels
    cross2D = centers(i,1)*normals(i,2) - centers(i,2)*normals(i,1);
    Cm = Cm - CP_vect(i) * lengths(i) / c^2 * cross2D;
end

end

