function [CP_vect] = CP_fun(U_inf,solution,geo,NPanels)
% This function calculates the CP distribution of an airfoil with the
% following formula: 
% 
% CP = 1 - U(x,y)^2 / |U_inf|^2

addpath mat_functions

% Compute elements from the panels 
[centers,~,tangent,extrema_1,extrema_2,~,~,L2G_TransfMatrix,G2L_TransfMatrix] = CreatePanels(geo);

% Compute CP 
CP_vect = zeros(NPanels,1);

for i = 1:NPanels

    local_center = centers(i, :)';
    U_panel = U_inf;

    for j = 1:NPanels
        local_extreme_1 = extrema_1(j, :)';
        local_extreme_2 = extrema_2(j, :)';
        local_L2G_TransfMatrix = squeeze(L2G_TransfMatrix(j, :, :));
        local_G2L_TransfMatrix = squeeze(G2L_TransfMatrix(j, :, :));
        U_panel = U_panel +...
            uSource(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix).*solution(j) + ...
            uVortex(local_center, local_extreme_1, local_extreme_2, local_L2G_TransfMatrix, local_G2L_TransfMatrix).*solution(end);
    end 
    U = dot(U_panel, tangent(i,:));

    % CP formula
    CP_vect(i) = 1 - (U / norm(U_inf))^2;
end

