%% WEISSINGER METHOD 
close all 
clear all
clc

%% Test Case 0: Wing of Cessna 172
% U_Inf_Mag = 62;     % [m/s]
% alpha = 3;          % [deg]
% beta = 0;           % [deg]
% U_Inf = [cosd(alpha)*cosd(beta) sind(beta) sind(alpha)] .* U_Inf_Mag;
% rho = 1.225;
% 
% config.NBodies = 1;
% 
% config.RootChord = [1.487];       % [m]
% config.DihedralAngle = [1.75];      % [deg]
% config.SweepAngle = [0];            % [deg]
% config.TaperRatio = [1];         
% config.Span = [10.89];            % [m]
% 
% config.LEPosition_X = [0];        % [m]
% config.LEPosition_Y = [0];
% config.LEPosition_Z = [0];        
% 
% config.RotationAngle_X = [0];       
% config.RotationAngle_Y = [0];   
% config.RotationAngle_Z = [0];
% 
% % Discretization options
% config.SemiSpanwiseDiscr = [20];
% config.ChordwiseDiscr = [20];
 

%% Test Case 1: Cessna 172
U_Inf_Mag = 62;     % [m/s]
alpha = 3;          % [deg]
beta = 0;           % [deg]
U_Inf = [cosd(alpha)*cosd(beta) sind(beta) sind(alpha)] .* U_Inf_Mag;
rho = 1.225;

config.NBodies = 2;

config.RootChord = [1.62, 0.7];       % [m]
config.DihedralAngle = [1.75, 0];      % [deg]
config.SweepAngle = [0, 0];            % [deg]
config.TaperRatio = [0.83, 1];         
config.Span = [11.0, 3.45];            % [m]

config.LEPosition_X = [0, 5.0];        % [m]
config.LEPosition_Y = [0, 0];
config.LEPosition_Z = [0, -1.1];        % tail is ~1.1 m below wing ref

config.RotationAngle_X = [0, 0];       
config.RotationAngle_Y = [0, -2];   % tail typically 0° incidence
config.RotationAngle_Z = [0, 0];

% Discretization options
config.SemiSpanwiseDiscr = [20,10];
config.ChordwiseDiscr = [20,10];

%% Test Case 2: Eurofigher Typhoon
% U_Inf_Mag = 62;     % [m/s]
% alpha = 3;          % [deg]
% beta = 0;           % [deg]
% U_Inf = [cosd(alpha)*cosd(beta) sind(beta) sind(alpha)] .* U_Inf_Mag;
% rho = 1.225;
% config.NBodies = 2;
% 
% % Wing (delta planform approximation)
% config.RootChord  = [8.31, 1.6];
% config.TaperRatio = [0.126, 0.4];
% config.Span       = [10.95, 3.0];
% config.SweepAngle = [46, 50];
% config.DihedralAngle = [0, -8];
% 
% config.LEPosition_X    = [0, -4.5];      % m
% config.LEPosition_Y    = [0, 0];
% config.LEPosition_Z    = [0, 0.5];       % m (tail below wing reference)
% 
% config.RotationAngle_X = [0, 0];
% config.RotationAngle_Y = [0, 2];         % tail incidence slightly lower
% config.RotationAngle_Z = [0, 0];
% 
% % Discretization options
% config.SemiSpanwiseDiscr = [20,20];
% config.ChordwiseDiscr = [20,20];

%% Preliminary computations

% Computing the span
config.SemiSpan = config.Span./2;
% Computing the surface
config.Surface = 2 * (config.SemiSpan .* config.RootChord .* ( 1 + config.TaperRatio ) ./ 2);
config.SurfaceProjected = config.Surface .* cosd(config.DihedralAngle);
% Computing the Tip chord
config.TipChord = config.RootChord .* config.TaperRatio;

% Compute MAC
config.MAC = (2/3) .* config.RootChord .* ( (1 + config.TaperRatio + config.TaperRatio.^2)./(1 + config.TaperRatio));

% Initialization of some variables --> to be consistent with literature
M = config.SemiSpanwiseDiscr;
N = config.ChordwiseDiscr;

% Dependent variables 
config.AspectRatio = config.Span.^2 ./ config.Surface;


%% Create the geometry structure

ControlPoints = cell(config.NBodies, 1);
InducedPoints = cell(config.NBodies, 1);
Normals = cell(config.NBodies, 1);
InfiniteVortices = cell(config.NBodies, 1);
Vortices = cell(config.NBodies, 1);
internalMesh = cell(config.NBodies, 1); 
WingExtremes = cell(config.NBodies, 1);


for iBody = 1:config.NBodies

    [ControlPoints{iBody}, InducedPoints{iBody}, Normals{iBody}, InfiniteVortices{iBody}, Vortices{iBody}, internalMesh{iBody}, WingExtremes{iBody}] = createStructure(config, iBody);

end
    

%% Matrices initialization

NPanelsTot = 2* config.SemiSpanwiseDiscr * config.ChordwiseDiscr';
matrixA = zeros(NPanelsTot, NPanelsTot);
knownTerm = zeros(NPanelsTot, 1);

%% Construction of the matrix

rowIndex = 0;

for iBody = 1:config.NBodies
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iBody)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iBody)
            
            % Update row index
            rowIndex = rowIndex + 1;
   
            columnIndex = 0;
            
            ControlPointHere = ControlPoints{iBody}{ChordPanel_i, SpanPanel_i}.Coords;
            LocalNormal = Normals{iBody}{ChordPanel_i, SpanPanel_i}.Coords;
            
            for jCorpo = 1:config.NBodies
                
                % Cycle on all of its chordwise panels
                for ChordPanel_j = 1:config.ChordwiseDiscr(jCorpo)
                    % Cycle on all of its spanwise panels
                    for SpanPanel_j = 1:2*config.SemiSpanwiseDiscr(jCorpo)
                        
                        % Update column index
                        columnIndex = columnIndex + 1;
                        
                        % Compute the influence induced by first
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.toInfty;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root.onWing;
                        U = vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by finite vortex
                        Extreme_1 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Root;
                        Extreme_2 = Vortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        % Compute the influence induced by second
                        % semi-infinite vortex
                        Extreme_1 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.onWing;
                        Extreme_2 = InfiniteVortices{jCorpo}{ChordPanel_j, SpanPanel_j}.Tip.toInfty;
                        U = U + vortexInfluence(ControlPointHere, Extreme_1, Extreme_2);

                        
                        matrixA(rowIndex, columnIndex) = dot(U, LocalNormal);
                       
                        
                    end
                end

            end
             
        end
    end
end

%% Known term construction 

rowIndex = 0;
for iBody = 1:config.NBodies
    
    % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iBody)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iBody)
            
            % Update row index
            rowIndex = rowIndex + 1;
  
            LocalNormal = Normals{iBody}{ChordPanel_i, SpanPanel_i}.Coords;
            
            knownTerm(rowIndex) = -dot(U_Inf, LocalNormal);
            
        end
    end
end

%% Solve the linear system

Solution = linsolve(matrixA, knownTerm);

Gamma = cell(config.NBodies, 1);
% Cells are useful if you're considering different bodies

rowIndex = 0;
for iBody = 1:config.NBodies
    
    Gamma{iBody} = zeros( config.ChordwiseDiscr(iBody), config.SemiSpanwiseDiscr(iBody)*2 );
    
     % Cycle on all of its chordwise panels
    for ChordPanel_i = 1:config.ChordwiseDiscr(iBody)
        % Cycle on all of its spanwise panels
        for SpanPanel_i = 1:2*config.SemiSpanwiseDiscr(iBody)
            
            % Update row index
            rowIndex = rowIndex + 1;
            
            Gamma{iBody}(ChordPanel_i, SpanPanel_i) = Solution(rowIndex);
        end
        
    end
    
end

% %% Visualization
% figure
% pcolor(Gamma{1}(:,:)); colorbar

%% Compute the 2D and 3D Lift
% 2D and Total Lift
[L_2D,L,delta_b] = Lift_fun(config,rho,U_Inf_Mag,Gamma,1);
% Compute CL 
config.CL = L/(1/2 * rho * U_Inf_Mag^2 * config.Surface(1));
fprintf('CL = %.3f\n',config.CL);

%% Compute 2D and 3D induced drag
% Compute control points for each 
CPoint = C4_ControlPoints(config,internalMesh);

% Compute induced drag 
[induced_alpha,D_ind_2D,D_ind] = Induced_Drag_fun(config,Normals,L_2D,CPoint,InfiniteVortices, U_Inf, Gamma,1);
config.CD_ind = D_ind/(1/2 * rho * U_Inf_Mag^2 * config.Surface(1));
fprintf('CDi = %.3f\n',config.CD_ind);

%% Visualization: 
% Plot of the aircraft
figure
hold on
axis equal
grid on
xlabel('X [m]')
ylabel('Y [m]')
zlabel('Z [m]')
title('Aircraft Geometry')

colors = lines(config.NBodies);

for iBody = 1:config.NBodies
    for j = 1:config.ChordwiseDiscr(iBody)
        for i = 1:2*config.SemiSpanwiseDiscr(iBody)

            % Panel corners
            LE_root = internalMesh{iBody}{j,i}.LERoot;
            LE_tip  = internalMesh{iBody}{j,i}.LEtip;
            TE_tip  = internalMesh{iBody}{j,i}.TEtip;
            TE_root = internalMesh{iBody}{j,i}.TERoot;

            X = [LE_root(1) LE_tip(1) TE_tip(1) TE_root(1)];
            Y = [LE_root(2) LE_tip(2) TE_tip(2) TE_root(2)];
            Z = [LE_root(3) LE_tip(3) TE_tip(3) TE_root(3)];

            patch(X, Y, Z, colors(iBody,:), ...
                  'FaceAlpha', 0.6, ...
                  'EdgeColor', 'k');
        end
    end
end

view(3)

% Color scale based on Gamma 
allGamma = [];
for iBody = 1:config.NBodies
    allGamma = [allGamma; Gamma{iBody}(:)];
end

cmin = min(allGamma);
cmax = max(allGamma);

colormap(jet)
clim([cmin cmax])
cb = colorbar;
cb.Label.String = '\Gamma [m^2/s]';

% Plot bound vortices
for iBody = 1:config.NBodies
    for j = 1:config.ChordwiseDiscr(iBody)
        for i = 1:2*config.SemiSpanwiseDiscr(iBody)

            V = Vortices{iBody}{j,i};

            x = [V.Root(1), V.Tip(1)];
            y = [V.Root(2), V.Tip(2)];
            z = [V.Root(3), V.Tip(3)];

            gamma_ji = Gamma{iBody}(j,i);

            % Normalize gamma for colormap
            t = (gamma_ji - cmin) / (cmax - cmin);
            t = max(0,min(1,t));
            color = interp1(linspace(0,1,64), jet(64), t);

            plot3(x,y,z,'LineWidth',2,'Color',color);
        end
    end
end

save config 
view(3)
