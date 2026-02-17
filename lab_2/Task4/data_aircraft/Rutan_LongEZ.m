function [config]=Rutan_LongEZ

config.NBodies = 2;        
config.RootChord = [2.37, 0.34];    
config.DihedralAngle = [1.7, 0]; 
config.SweepAngle1 = [47.5, 0];    

config.TaperRatio = [0.2194, 1];  
config.Span = [7.96, 3.50];                 

config.d=[1.53, 0]; 

% corda dove ho innesto del secondo tratto 
config.MediumChord=[1.05, 0];
 
config.SweepAngle2 = [22.7, 0]; 

config.LEPosition_X = [0.70, 0];   
config.LEPosition_Z = [0, 0];   
config.LEPosition_Y = [0, 0];   

config.RotationAngle_Y = input(['enter the value of the Angle of attack:\n' ...
    '(vector->[wing,tail], tipicaly on this aircraft a_tail=a_wing+1)):\n']);

% Discretization options
config.SemiSpanwiseDiscr = input('Choose number of pannel for the semi-Span discretization (vector->[wing,tail]):\n');
config.ChordwiseDiscr = input('Choose number of pannel for the Chordwise discretization (vector->[wing,tail]):\n');

%% Preliminary computations

% Computing the semi span
config.SemiSpan = config.Span(1:end)./2;

% Computing the surface
config.Surface = 2 * (config.SemiSpan(1:end) .* config.RootChord(1:end) .* ( 1 + config.TaperRatio(1:end) ) ./ 2); % x

% Computing the Tip chord
config.TipChord = config.RootChord(1:end) .* config.TaperRatio(1:end); 