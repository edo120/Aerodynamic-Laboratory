function [v] = biot_sa(Center,Extreme1,Extreme2) 

% This function computes the induced contribution using the integral Biot–Savart 
% formulation over a vortex segment

rc1 = Center - Extreme1;
r21 = Extreme2 - Extreme1;
rc2 = Center - Extreme2;

r1_norm=rc1/norm(rc1);
r2_norm=rc2/norm(rc2);

toll = 1e-8;
CP_norm = norm(cross(rc1,rc2));

if CP_norm < toll
    CP_norm = toll;
end

v = (1/(4*pi)) * dot(r21,(r1_norm-r2_norm)) * cross(rc1,rc2)'/((CP_norm)^2);


