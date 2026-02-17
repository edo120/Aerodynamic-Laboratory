function [polarData] = createData(airfoildata,NPanels,Re,alpha,Ncrit)
% This function uses xfoil to create the polar and extract 
% the coordinates of the points 


% creating/rewriting working directory
workDir = fullfile(pwd,'xfoil_workdir_task2');
if exist(workDir,'dir')
    rmdir(workDir,'s');
end
mkdir(workDir);


% initialising
n_Re = numel(Re);
n_Ncrit = numel(Ncrit);
a = numel(alpha);

polarData = struct([]);

% cycles
for i = 1 : n_Re
    for j = 1 : n_Ncrit

        Re_i = Re(i);
        Ncrit_j = Ncrit(j);

        % Matlab generates a script with the polar and the instruction for xfoil
        xfoilPrompt = fullfile(workDir,'xfoil.prompt');
        polarFile = sprintf('polar_Re%d_N%d.dat',Re_i,Ncrit_j);
        boundary_layer = string(a);

        fileID = fopen(xfoilPrompt,'w');
        fprintf(fileID, 'load %s\n',airfoildata);
        fprintf(fileID,'ppar\n');
        fprintf(fileID,'n %d\n\n\n',NPanels);
        fprintf(fileID,'mdes\n');
        fprintf(fileID,'filt\n');
        fprintf(fileID,'exec \n\n');
    
        fprintf(fileID,'pane\n');
        fprintf(fileID,'oper\n');
        fprintf(fileID,'iter 500\n');
        fprintf(fileID,'visc %e\n',Re_i);
        fprintf(fileID,'vpar\n');
        fprintf(fileID,'N %d\n\n',Ncrit_j);
        fprintf(fileID,'pacc\n');
        fprintf(fileID,'%s\n\n',polarFile);

        for k = 1 : a
            fprintf(fileID,'alfa %f\n',alpha(k));
            boundary_layer(k) = sprintf('bl_Re%d_N%d_%f',Re_i,Ncrit_j,alpha(k));
            fprintf(fileID,'dump\n');
            fprintf(fileID,'%s\n',boundary_layer(k));
        end

        fprintf(fileID,'pacc\n\n');
        fprintf(fileID,'quit\n');

        fclose(fileID);


        % running xfoil
        system(sprintf('xfoil < %s > /dev/null',xfoilPrompt));
                

        % reading results
        data = readPolarFile(polarFile);

        
        % saving
        polarData(i,j).Re = Re_i;
        polarData(i,j).Ncrit = Ncrit_j;
        polarData(i,j).data = data;


        % moving the file in xfoil_workdir
        movefile(polarFile,'xfoil_workdir_task2/');
        for k = 1 : a
            bl = read_bl(boundary_layer(k));
            polarData(i,j).data_bl.alpha(k,1) = bl;
            movefile(boundary_layer(k),'xfoil_workdir_task2/');
        end

        if exist(':00.bl','file')
            movefile(':00.bl','xfoil_workdir_task2/');
        end
        
    
    end
end
end