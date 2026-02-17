function corpo = importXfoilProfile(filename) 
% This function reads Xfoil file and return data as a struct
%
%   corpo.x
%   corpo.y

% Verify file opening
if ~isfile(filename)
    error('Polar file not found: %s', filename);
end

% Reading data skipping Xfoil's Header
data = readmatrix(filename, ...
    'FileType','text', ...
    'NumHeaderLines',1);

% Check
if size(data,2) < 2
    error('Invalid polar format (insufficients columns).');
end

corpo.x = data(:,1);
corpo.y = data(:,2);

end

