function bl = read_bl(filename)
% This function reads Xfoil boundary layer file and return data as a struct
%
% bl.s      
% bl.x     
% bl.y      
% bl.Ue_Vinf 
% bl.Dstar  
% bl.Theta   
% bl.Cf    
% bl.H     
% bl.H_star  
% bl.p       
% bl.m       
% bl.k       

% Verify file opening
if ~isfile(filename)
    error('Polar file not found: %s', filename);
end

% Reading data skipping Xfoil's Header
data = readmatrix(filename, ...
    'FileType','text', ...
    'NumHeaderLines',1);
%data = data(:,2:end);

% Check
if size(data,2) < 12
    error('Invalid polar format (insufficients columns).');
end

% Struct
bl.s       = data(:,1);
bl.x       = data(:,2);
bl.y       = data(:,3);
bl.Ue_Vinf = data(:,4);
bl.Dstar   = data(:,5);
bl.Theta   = data(:,6);
bl.Cf      = data(:,7);
bl.H       = data(:,8);
bl.H_star  = data(:,9);
bl.p       = data(:,10);
bl.m       = data(:,11);
bl.k       = data(:,12);

end