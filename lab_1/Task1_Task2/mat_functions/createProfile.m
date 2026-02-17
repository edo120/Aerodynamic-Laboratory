function [x,y,polarData] = createProfile(naca,NPanels,Chord,alpha)
% This function uses xfoil to create the profile and extract 
% the coordinates of the points of the profile
% for the Hess-Smith algorithm

% creating/rewriting working directory
workDir = fullfile(pwd,'xfoil_workdir_task1');
if exist(workDir,'dir')
    rmdir(workDir,'s');
end
mkdir(workDir);

a = numel(alpha);

% initialising string names
xfoilPrompt = fullfile(workDir,'xfoil_prompt_task1.txt');
filename = sprintf('NACA_%s.dat',naca);
polarFile = sprintf('polar_NACA_%s.dat',naca);
Cp_distribution = string(a);


% Matlab generates a script with the instruction for xfoil
fileID = fopen(xfoilPrompt,'w');
fprintf(fileID, 'naca %s\n',naca);
fprintf(fileID,'pane\n');
fprintf(fileID,'gdes\n');
fprintf(fileID,'tgap 0 0 \n');
fprintf(fileID,'exec \n\n');
fprintf(fileID,'ppar\n');
fprintf(fileID, 'n %d\n\n\n',NPanels+1); 
fprintf(fileID, 'save %s\n',filename);

fprintf(fileID,'pane\n');
fprintf(fileID,'xycm 0.0 0.0\n');  % change the pole of Cm: AC-->LE
fprintf(fileID,'oper\n');
fprintf(fileID,'iter 500\n');
fprintf(fileID,'pacc\n');
fprintf(fileID,'%s\n\n',polarFile);
for k = 1 : a
    fprintf(fileID,'alfa %f\n',alpha(k));
    Cp_distribution(k) = sprintf('Cp_distribution_NACA_%s_%f',naca,alpha(k));
    fprintf(fileID,'cpwr\n');
    fprintf(fileID,'%s\n',Cp_distribution(k));
end
fprintf(fileID,'pacc\n\n');
fprintf(fileID,'quit \n');
fclose(fileID);


% running xfoil
system(sprintf('xfoil < %s > /dev/null',xfoilPrompt));


% reading results      
Corpo = importXfoilProfile(filename);
polarData = readPolarFile(polarFile);

% Invert vectors
x = flipud(Corpo.x);
y = flipud(Corpo.y);

x = x.*Chord;
y = y.*Chord;


% moving the file in xfoil_workdir
movefile(filename,'xfoil_workdir_task1/');
movefile(polarFile,'xfoil_workdir_task1/');
for k = 1 : a
    Cp = importXfoilProfile(Cp_distribution(k));
    polarData.data_Cp.alpha(k,1) = Cp;
    movefile(Cp_distribution(k),'xfoil_workdir_task1/');
end
if exist(':00.bl','file')
    movefile(':00.bl','xfoil_workdir_task1/');
end

end