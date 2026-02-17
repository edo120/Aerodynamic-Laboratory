function [config]=Concorde_data

config.NBodies = 1;

config.RootChord  = 27.66;
config.TaperRatio = 0.17;
config.Span= 25.56;
config.SweepAngle1 = 0;
config.DihedralAngle = -1.7;

config.LEPosition_X    = 0;   
config.LEPosition_Y    = 0;
config.LEPosition_Z    = 0;   

config.RotationAngle_X = 0;
config.RotationAngle_Y = input('enter the value of the Angle of attack:\n');
config.RotationAngle_Z = 0;

% Discretization options
config.SemiSpanwiseDiscr = input('Choose number of pannel for the semi-Span discretization:\n');
config.ChordwiseDiscr = input('Choose number of pannel for the Chordwise discretization:\n');

% Preliminary computations

% Computing the Tip chord
config.TipChord = config.RootChord(1:end) .* config.TaperRatio(1:end); 

% per avere ala a delta dritta dietro 
config.SweepAngle1(1)=atand((config.RootChord(1)*(1-config.TaperRatio(1)))/(config.Span(1)/2));

% Computing the semi span
config.SemiSpan = config.Span(1:end)./2;

% Computing the surface
config.Surface = 2 * (config.SemiSpan(1:end) .* config.RootChord(1:end) .* ( 1 + config.TaperRatio(1:end) ) ./ 2); % x

config.d=0;

config.SweepAngle2 = 0; 

config.MediumChord = (config.RootChord(1)+config.TipChord(1))/2 ;