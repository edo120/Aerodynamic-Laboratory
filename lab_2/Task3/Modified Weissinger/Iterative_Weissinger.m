clear 
close all 
clc

tic

alfa=[-2:2];
data.plot_flag=0;

data.U_mod=62;
data.alpha=0;
data.beta=0; 
data.rho=1.225;

data.NBodies = 2;
data.SemiSpanwiseDiscr = [20,10];
data.ChordwiseDiscr = [20,10];

data.RootChord = [1.487, 1.3];   
data.DihedralAngle = [1.75, 0]; 
data.SweepAngle1 = [0, 0];   

data.TaperRatio = [1, 0.615];  
data.Span = [10.89, 3.5];         

data.LEPosition_X = [0, 4.5];   
data.LEPosition_Z = [0, -0.5];   
data.LEPosition_Y = [0, 0];   

data.RotationAngle_Y = [3,1];

data.d=[0, 0];
data.MediumChord=[0, 0];
data.SweepAngle2 = [0, 0]; 

SemiSpan = data.Span(1:end)./2;

data.SweepAngle1(2) = [(atand((1-data.TaperRatio(2))*data.RootChord(2)/(SemiSpan(2)-data.d(2))))/2];

[~, ~, ~, ~, Cl, Cd]=weissinger_fun(alfa, data);

toc
confront.Cl=Cl;
confront.Cd=Cd;
save Wessinger_total2 confront

%%

tic
alfa=[-3:2];
data.plot_flag=0;

data.U_mod=62;
data.alpha=0;
data.beta=0; 
data.rho=1.225;

data.NBodies = 1;
data.SemiSpanwiseDiscr = [20];
data.ChordwiseDiscr = [20];

data.RootChord = [1.487, 1.3];   
data.DihedralAngle = [1.75, 0]; 
data.SweepAngle1 = [0, 0];   

data.TaperRatio = [1, 0.615];  
data.Span = [10.89, 3.5];         

data.LEPosition_X = [0, 4.5];   
data.LEPosition_Z = [0, -0.5];   
data.LEPosition_Y = [0, 0];   

data.RotationAngle_Y = [3,1];

data.d=[0, 0];
data.MediumChord=[0, 0];
data.SweepAngle2 = [0, 0]; 

SemiSpan = data.Span(1:end)./2;

data.SweepAngle1(2) = [(atand((1-data.TaperRatio(2))*data.RootChord(2)/(SemiSpan(2)-data.d(2))))/2];

[~, ~, ~, ~, Cl, Cd]=weissinger_fun(alfa, data);
toc

confront.Cl=Cl;
confront.Cd=Cd;
save Wessinger_single2 confront