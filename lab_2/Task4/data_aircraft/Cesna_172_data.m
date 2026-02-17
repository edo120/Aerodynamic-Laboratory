function [config]=Cesna_172_data

config.NBodies = input(['Enter number of body:\n' ...
    '1 for just the wing\n' '2 for the tow-surface configuration\n']); 

config.RootChord = [1.62, 1.3];   
config.DihedralAngle = [1.44, 0];  
config.SweepAngle1 = [0, 0];    

% corda dove ho innesto del secondo tratto 
config.MediumChord=[config.RootChord(1), 0];

config.TaperRatio = [0.687, 0.615];  
config.Span = [10.92, 3.45];        

config.LEPosition_X = [0, 4.5];    
config.LEPosition_Z = [0, -0.5];   
config.LEPosition_Y = [0, 0];   

config.RotationAngle_Y = input(['enter the value of the Angle of attack:\n' ...
    '(vector->[wing,tail], tipicaly on this aircraft a_tail=a_wing-2)\n']);

% Discretization options
config.SemiSpanwiseDiscr = input('Choose number of pannel for the semi-Span discretization (for 2 surface vector->[wing,tail]):\n');
config.ChordwiseDiscr = input('Choose number of pannel for the Chordwise discretization (for 2 surface vector->[wing,tail]):\n');

%% Preliminary computations

% Computing the semi span
config.SemiSpan = config.Span(1:end)./2;

% Computing the surface
config.Surface = 2 * (config.SemiSpan(1:end) .* config.RootChord(1:end) .* ( 1 + config.TaperRatio(1:end) ) ./ 2); % x

% Computing the Tip chord
config.TipChord = config.RootChord(1:end) .* config.TaperRatio(1:end); % x

% parametro utile in un futuro 
config.d=[0.273*12, 0]; % lunghezza tratto retilineo nella semiala 


 
config.SweepAngle2 = [0, 0]; 

config.SweepAngle1(2) = (atand((1-config.TaperRatio(2))*config.RootChord(2)/(config.SemiSpan(2)-config.d(2))))/2;  

 
