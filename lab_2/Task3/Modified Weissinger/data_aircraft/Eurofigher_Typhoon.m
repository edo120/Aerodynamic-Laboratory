function [config]=Eurofigher_Typhoon

config.NBodies = 2;

config.RootChord  = [8.31, 1.6];
config.TaperRatio = [0.126, 0.4];
config.Span = [10.95, 3.0];
config.SweepAngle1 = [0, 50];
config.DihedralAngle = [0, -8];

config.LEPosition_X    = [0, -4.5];      
config.LEPosition_Y    = [0, 0];
config.LEPosition_Z    = [0, 0.5];     

config.RotationAngle_Y = input('enter the value of the Angle of attack (vector->[wing,tail]):\n');

config.SemiSpanwiseDiscr = input('Choose number of pannel for the semi-Span discretization (vector->[wing,tail]):\n');
config.ChordwiseDiscr = input('Choose number of pannel for the Chordwise discretization (vector->[wing,tail]):\n');

%% Preliminary computations

% Computing the Tip chord
config.TipChord = config.RootChord(1:end) .* config.TaperRatio(1:end); 

% per avere ala a delta dritta dietro 
config.SweepAngle1(1)=atand((config.RootChord(1)*(1-config.TaperRatio(1)))/(config.Span(1)/2));

% Computing the semi span
config.SemiSpan = config.Span(1:end)./2;

% Computing the surface
config.Surface = 2 * (config.SemiSpan(1:end) .* config.RootChord(1:end) .* ( 1 + config.TaperRatio(1:end) ) ./ 2); % x

config.d=[0, 0];  

config.MediumChord=[(config.RootChord(1)+config.TipChord(1))/2, (config.RootChord(2)+config.TipChord(2))/2];

config.SweepAngle2 = [0, 0];  