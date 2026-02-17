function [v] = biot_sa2(Center,Extreme1,Eextreme2, Tangent, Extreme_inf)

% This function computes the induced contribution using the integral Biot–Savart 
% formulation over a vortex segment, sending one end of the vortex to infinity, 
% if needed, with the same inclination as the single panel.

if nargin > 4
    if Extreme_inf == Eextreme2

       Extreme1 = Eextreme2;
       Eextreme2 = 1e3*Tangent + Eextreme2; 

    elseif Extreme_inf == Extreme1
       
       Eextreme2 = Extreme1;
       Extreme1 = 1e3*Tangent + Extreme1;
       
    end
else
end


rc1 = Center - Extreme1;
r21 = Eextreme2 - Extreme1;
rc2 = Center - Eextreme2;    

r1_norm=rc1/norm(rc1);
r2_norm=rc2/norm(rc2);

toll = 1e-8;
CP_norm= norm(cross(rc1,rc2));

if CP_norm < toll
    CP_norm = toll;
end

v = (1/(4*pi)) * dot(r21,(r1_norm-r2_norm)) * cross(rc1,rc2)'/((CP_norm)^2);


