function [Cl] = Cl_fun(U_inf,geo,CP_vect,NPanels)
% This function calculates the Cl (lower case l --> 2D airfoil) of an
% airfoil

addpath mat_functions

% Compute elements from the panels 
[~,normals,~,~,~,~,lengths,~,~] = CreatePanels(geo);

% Compute Cl 

Cl = 0;
c = max(geo.x) - min(geo.x);                  % Chord length
tangent_U = U_inf/norm(U_inf);                % Tangential vector of U_inf
normal_U = [-tangent_U(2),tangent_U(1)];      % Normal vector to U_inf

for i = 1:NPanels
    Cl = Cl - (CP_vect(i)*lengths(i)/c) * dot(normals(i,:),normal_U); 
end

