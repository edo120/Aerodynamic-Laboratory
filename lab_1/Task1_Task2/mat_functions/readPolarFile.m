function polar = readPolarFile(filename)
% This function reads Xfoil polar file and return data as a struct
%
%   polar.alpha   [deg]
%   polar.CL
%   polar.CD
%   polar.CDp
%   polar.CM
%   polar.Top_Xtr
%   polar.Bot_Xtr

% Verify file opening
if ~isfile(filename)
    error('Polar file not found: %s', filename);
end

% Reading data skipping Xfoil's Header
data = readmatrix(filename, ...
    'FileType','text', ...
    'NumHeaderLines',12);

% Check
if size(data,2) < 7
    error('Invalid polar format (insufficients columns).');
end

% Struct
polar.alpha   = data(:,1);
polar.Cl      = data(:,2);
polar.Cd      = data(:,3);
polar.Cdp     = data(:,4);
polar.Cm      = data(:,5);
polar.Top_Xtr = data(:,6);
polar.Bot_Xtr = data(:,7);


end
