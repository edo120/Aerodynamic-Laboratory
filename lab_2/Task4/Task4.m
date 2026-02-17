clear 
close all
clc

% If you want to perform a quicker check, you can flag the multiple-angles section; 
% the data are already saved, and by running the plotting part they will be loaded automatically

data_cessna.U_mod=62.8;
data_cessna.alpha=0;
data_cessna.beta=0; 
data_cessna.rho=1.225;

alfa=[0];
data_cessna.plot_flag=1;
data_cessna.NBodies = 2;         
data_cessna.RootChord = [1.63, 1.3];   
data_cessna.DihedralAngle = [1.44, 0]; 
data_cessna.SweepAngle1 = [0, 0];   

data_cessna.TaperRatio = [0.687, 0.615];  
data_cessna.Span = [10.92, 3.45];         

data_cessna.LEPosition_X = [0, 4.5];   
data_cessna.LEPosition_Z = [0, -0.5];   
data_cessna.LEPosition_Y = [0, 0];   

data_cessna.RotationAngle_Y = [3,-1];

data_cessna.SemiSpanwiseDiscr = [20,10];
data_cessna.ChordwiseDiscr = [20,10];
 
data_cessna.d=[0.273*12, 0];
data_cessna.MediumChord=[data_cessna.RootChord(1), 0];
data_cessna.SweepAngle2 = [0, 0]; 

SemiSpan = data_cessna.Span(1:end)./2;

data_cessna.SweepAngle1(2) = (atand((1-data_cessna.TaperRatio(2))*data_cessna.RootChord(2)/(SemiSpan(2)-data_cessna.d(2))))/2;

[output_flux_cessna, wing_data_cessna, sol_cessna, imp_points_cessna, Cl_cessna, Cd_cessna ]=weissinger_fun(alfa, data_cessna);


% rutan con angolo di incidenza dell'ala principale 3 gradi 

alfa=[0];

data_rutan.U_mod=62.8;
data_rutan.alpha=0;
data_rutan.beta=0; 
data_rutan.rho=1.225;

data_rutan.plot_flag=1;

data_rutan.NBodies = 2;        
data_rutan.RootChord = [2.37, 0.34];    
data_rutan.DihedralAngle = [1.7, 0]; 
data_rutan.SweepAngle1 = [47.5, 0];    

data_rutan.TaperRatio = [0.2194, 1];  
data_rutan.Span = [7.96, 3.50];                 

data_rutan.d=[1.53, 0]; 
data_rutan.MediumChord=[1.05, 0];
 
data_rutan.SweepAngle2 = [22.7, 0]; 
data_rutan.LEPosition_X = [0.70, 0];   
data_rutan.LEPosition_Z = [0, 0];   
data_rutan.LEPosition_Y = [0, 0];   

data_rutan.RotationAngle_Y = [3,4];

data_rutan.SemiSpanwiseDiscr = [20,10];
data_rutan.ChordwiseDiscr = [20,10];

[output_flux_rutan, wing_data_rutan, sol_rutan, imp_points_rutan,Cl_rutan, Cd_rutan ]=weissinger_fun(alfa, data_rutan);

%%
figure 
xx2 = linspace(-data_cessna.Span(1)/2,data_cessna.Span(1)/2,data_cessna.SemiSpanwiseDiscr(1)*2);
plot(xx2, sol_cessna.ind_alpha{1,1}, 'o-','LineWidth',1.5)
hold on
plot(xx2, sol_rutan.ind_alpha{1,1},'o-','LineWidth',1.5)
grid on 
xlabel('Spanwise position [m]', FontSize= 12)
ylabel('\alpha [deg]', "Rotation", 0, FontSize=12)
legend({'wing cessna 172','wing rutan-log ez'},'Location', 'north')
title('Comparison of Wing Induced Angles for Selected Aircraft Models', FontSize=12)

%% multi angles
% 
% data_cessna.plot_flag=0;
% alfa=[-3:3];
% [~, ~, ~, ~,Cl_cessna_multi, Cd_cessna_multi ]=weissinger_fun(alfa, data_cessna);
% 
% 
% 
% data_rutan.plot_flag=0;
% [~, ~, ~, ~,Cl_rutan_multi, Cd_rutan_multi ]=weissinger_fun(alfa, data_rutan);
% 
% %%
% 
% confront.Cl_cessna=Cl_cessna_multi;
% confront. Cd_cessna= Cd_cessna_multi;
% confront.Cl_rutan=Cl_rutan_multi;
% confront. Cd_rutan=Cd_rutan_multi;
% 
% save Task4_multi_sol confront

%%

load Task4_multi_sol.mat
alfa=[-3:3];

Cl_cessna_multi=confront.Cl_cessna;
Cd_cessna_multi=confront. Cd_cessna;
Cl_rutan_multi=confront.Cl_rutan;
Cd_rutan_multi=confront. Cd_rutan;

figure
plot(Cd_rutan_multi, Cl_rutan_multi, 'o-','LineWidth',1.5)
hold on 
grid on
plot(Cd_cessna_multi, Cl_cessna_multi, 'o-','LineWidth',1.5)
xlabel('CD', FontSize=12)
ylabel('CL', "Rotation", 0, FontSize=12)
title('Polar', FontSize=18)
legend('rutan-log ez', 'cessna 172', 'Location', 'northwest')


figure 
plot(alfa, Cl_rutan_multi, 'o-','LineWidth',1.5)
hold on 
grid on
plot(alfa, Cl_cessna_multi, 'o-','LineWidth',1.5)
xlabel('\alpha [deg]', FontSize=12)
ylabel('CL', "Rotation", 0, FontSize=12)
title('CL/\alpha curve', FontSize=18)
legend('rutan-log ez', 'cessna 172', 'Location', 'northwest')

