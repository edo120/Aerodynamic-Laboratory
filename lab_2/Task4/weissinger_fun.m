function [output_flux, wing_data, sol, imp_points, Cl_zita, Cd_zita ]=weissinger_fun(zita, data)

addpath data_aircraft
addpath  weissinger_fun

if nargin<2
s=input(['choose your aircraft wing from a selection, semi-manual or manual mod\n' 'enter\n' '1 for data base setting\n' ...
    '2 for semi-manual setting\n']);

    if s ==1
        aircraft=input(['1-->Cesna 172\n' '2-->Concorde\n' '3-->Eurofigher Typhoon\n' '4-->Rutan Long-EZ\n']);
    
        if aircraft==1
            [config]=Cesna_172_data;
        end
        if aircraft==2
            [config]=Concorde_data;
        end
        if aircraft==3
            [config]=Eurofigher_Typhoon;
        end
        if aircraft==4
            [config]=Rutan_LongEZ;
        end
    
    elseif s==2
    % inserimento dei dati 
    
    input(['PRESS ENTER AFTER READING TO START\n'...
        'You have selected the program’s semi-manual mode; therefore, you must enter the characteristic data of your wing.\n'...
        'If the aircraft includes multiple surfaces, enter the data as vectors (e.g., [wing, tail]).\n' ...
        'The sweep angles are measured from the horizontal.'])
    
    config.NBodies= input('enter the number of bodies that make up the configuration:\n');
    config.RootChord  = input('enter the value of the Root chord:\n');
    config.TaperRatio = input('enter the value of the Taper ratio:\n');
    config.Span = input('enter the value of the Span:\n');
    config.SweepAngle1 = input('enter the value of the Sweep Angle 1:\n');
    config.DihedralAngle = input('enter the value of the Dihedral Angle:\n');
    
    % parametro utile in un futuro 
    config.d= input('enter the length of the rectangular section along the semi-span.:\n');
    
    config.MediumChord=input(['enter the value of the chord in the point d:\n' ...
        '!! Enter only if the body (wing or tail) has a double geometry (rectangular–trapezoidal), otherwise, enter 0\n']);
    
    config.SweepAngle2 = input(['enter the value of the Sweep Angle 2:\n' ...
        '!! Enter only if the body (wing or tail) has a double geometry (rectangular–trapezoidal), otherwise, enter 0\n']);
    
    config.LEPosition_X    = input('enter the position of the body in x direction:\n');
    config.LEPosition_Y    = input('enter the position of the body in y direction (probably [0,0]):\n');
    config.LEPosition_Z    = input('enter the position of the body in z direction:\n');
    
    config.RotationAngle_Y = input('enter the value of the wing incidence:\n');
    
    config.SemiSpanwiseDiscr = input('Choose number of pannel for the semi-Span discretization:\n');
    config.ChordwiseDiscr = input('Choose number of pannel for the Chordwise discretization:\n');
    
    
    % Computing the Tip chord
    config.TipChord = config.RootChord(1:end) .* config.TaperRatio(1:end); 
    
    % Computing the semi span
    config.SemiSpan = config.Span(1:end)./2;
    
    % Computing the surface
    config.Surface = 2 * (config.SemiSpan(1:end) .* config.RootChord(1:end) .* ( 1 + config.TaperRatio(1:end) ) ./ 2); % x 
    end

data_flux.U_mod = input('flux velocity module:\n');
data_flux.alpha = input('enter incidence flux angle:\n');
data_flux.beta = input('enter sideslip angle:\n'); % sideslip angle 
data_flux.rho = 1.225;


plot_flag=input(['Do you want the plots as output?\n' '1 for yes\n' '0 for no\n']);
end


if nargin==2
    config.NBodies=data.NBodies;
    config.RootChord  = data.RootChord;
    config.TaperRatio = data.TaperRatio;
    config.Span = data.Span;
    config.SweepAngle1 = data.SweepAngle1;
    config.DihedralAngle = data.DihedralAngle;

    config.d=data.d;
    config.MediumChord=data.MediumChord;
    config.SweepAngle2 = data.SweepAngle2;

    config.LEPosition_X    =  data.LEPosition_X;
    config.LEPosition_Y    = data.LEPosition_Y;
    config.LEPosition_Z    = data.LEPosition_Z;
    
    config.RotationAngle_Y =  data.RotationAngle_Y;
    
    config.SemiSpanwiseDiscr =  data.SemiSpanwiseDiscr;
    config.ChordwiseDiscr =  data.ChordwiseDiscr;
      
    % Computing the Tip chord
    config.TipChord = config.RootChord(1:end) .* config.TaperRatio(1:end); 
    
    % Computing the semi span
    config.SemiSpan = config.Span(1:end)./2;
    
    % Computing the surface
    config.Surface = 2 * (config.SemiSpan(1:end) .* config.RootChord(1:end) .* ( 1 + config.TaperRatio(1:end) ) ./ 2); % x

    data_flux.U_mod = data.U_mod;
    data_flux.alpha = data.alpha;
    data_flux.beta = data.beta; 
    data_flux.rho = data.rho;

    plot_flag=data.plot_flag;

end

%%
    Cl_zita=zeros(length(zita), 1);
    Cd_zita=zeros(length(zita), 1);

    AoA=config.RotationAngle_Y;

for f=1:length(zita)
    
    config.RotationAngle_Y=AoA+zita(f);
    %%
    structure_data = struct('choord_panels',config.ChordwiseDiscr,'semiSpan_panels',config.SemiSpanwiseDiscr,...
        'LEposition_x',config.LEPosition_X,'LEposition_z',config.LEPosition_Z, ...
        'RotationA_y',config.RotationAngle_Y);
    
    wing_data = struct('diedro',config.DihedralAngle,'sweep1',config.SweepAngle1,...
        'sweep2',config.SweepAngle2,'b', config.Span,'Sup',config.Surface,...
        'cr',config.RootChord, 'ct',config.TipChord,'d',config.d, 'cmw', config.MediumChord);
    
    %% 
    
    % creo le directory dove salverò i dati 
    
    imp_points = struct('E1', {{}}, 'E2', {{}}, 'Center', {{}},'Normal',{{}},'W',{{}},'Choord',{{}});
    plot_data = struct('X', {{}}, 'Y', {{}}, 'Z', {{}},'Ex',{{}},'Ey',{{}},'Ez',{{}},'Cx',{{}},'Cy',{{}},'Cz',{{}});
    
    % qua sotto creo paneli
    
    for z=1:config.NBodies 
        [imp_points,plot_data]=create_Geometry(structure_data,wing_data,imp_points,plot_data,z);
    end
    
    N_pannels_x = structure_data.choord_panels;        
    N_pannels_y = 2 * structure_data.semiSpan_panels; 
    N_pannels_bodies = N_pannels_x .* N_pannels_y;   

    %% create output flux
    
    alfa=data_flux.alpha;
    beta=data_flux.beta;
    U_inf_mod=data_flux.U_mod;
    
    alfa = deg2rad(alfa * ones(1,config.NBodies));
    beta = deg2rad(beta * ones(1,config.NBodies));
    
    for z = 1:config.NBodies
    
        U_inf = [cos(alfa(z))*cos(beta(z)), sin(beta(z))*cos(alfa(z)), sin(alfa(z))] * U_inf_mod;
        output_flux.U_inf{z} = U_inf;
        output_flux.alpha{z} = alfa(z);
        output_flux.Rxa{z} = [cos(alfa(z)) 0 sin(alfa(z)); 0 1 0; -sin(alfa(z)) 0 cos(alfa(z))]*...
            [cos(beta(z)) -sin(beta(z)) 0; sin(beta(z)) cos(beta(z)) 0; 0 0 1 ];
    
    end
    
    output_flux.rho = data_flux.rho;
    
    %%                        
    
    for z = 1 : config.NBodies         
    
        N_pannels_body1 = N_pannels_bodies(z);     
    
        for k = 1: config.NBodies    
    
            N_pannels_body2 = N_pannels_bodies(k);           
            A = zeros(N_pannels_body1,N_pannels_body2);     
            KnowTerms = zeros(N_pannels_body1,1);           
            v_ind_mat = zeros(N_pannels_body1,N_pannels_body2);
    
            for l = 1 : N_pannels_body1 
    
                Control_point = imp_points.Center{z}(l,:);         
                Normal = imp_points.Normal{z}(l,:);          
                KnowTerms(l) = -dot(output_flux.U_inf{z},Normal);  

                for r = 1 : N_pannels_body2                
                    
                    T=imp_points.Tangent{k}(r, :);     % estraggo la tantente relativa al panello da cui partirà il vortice
                    Extreme_1 = imp_points.E1{k}(r,:);        % estraggo E1 a righe cordinate 
                    Extreme_2 = imp_points.E2{k}(r,:);        % estraggo E2 a righe cordinate

                    v_Einf1 = biot_sa2(Control_point,Extreme_1,Extreme_2,T, Extreme_1);       
                    v_E1E2 = biot_sa(Control_point,Extreme_1,Extreme_2);    
                    v_Einf2 = biot_sa2(Control_point,Extreme_1,Extreme_2,T, Extreme_2);

                    v_ind=  v_Einf1 + v_Einf2;
                    v_ind_mat(l,r) = (v_ind(3)*imp_points.Normal{k}(r, 3));    
    
                    V_gamma = v_E1E2 + v_Einf1 + v_Einf2;  % nel file è q_jk

                    A(l,r) = dot(V_gamma,Normal);

                end
                
                matrix_single(z).A_scomp{k} = A;
                matrix_single(z).v_ind{k} = v_ind_mat;
               
            end
        end
        
        matrix_single(z).KT{1} = KnowTerms;
    end
    
   
    A=[];
    KnowTerms=[];
    v_ind_mat = [];

    for z = 1: config.NBodies     
        A_scomp = cell2mat(matrix_single(z).A_scomp);
        KT = cell2mat(matrix_single(z).KT);
        v_ind = cell2mat(matrix_single(z).v_ind);
        A = [A;A_scomp];
        KnowTerms = [KnowTerms;KT];
        v_ind_mat = [v_ind_mat;v_ind]; 
    end
    
    GAMMA = A\KnowTerms;                
    V_INFP = v_ind_mat* GAMMA;     
    
    clear A KnowTerms v_ind_mat 
    
    %% POST PROCESSING
    
    i_p =0;
    
    for z = 1 : config.NBodies
       
        N_pannels_body1 = N_pannels_bodies(z);   
        gamma= GAMMA(i_p + 1 : i_p + N_pannels_body1);
    
        v_infp = V_INFP(i_p + 1 : i_p + N_pannels_body1);
    
        s = N_pannels_y(z);
        gamma_sp = zeros(N_pannels_x(z),N_pannels_y(z));

        for k = 1:s
            gamma_sp(:,k) = gamma(k:s:N_pannels_body1);
        end

        sol.gamma_sp{z} = gamma_sp;
    
        Aerodynamic_forces = zeros(N_pannels_body1,3);
        
        
        for i = 1: N_pannels_body1
            U_vinf=[data_flux.U_mod, 0, v_infp(i)]; 
            Aerodynamic_forces(i,:) = -output_flux.rho*cross([0,1,0]*gamma(i),U_vinf);
        end
        
        L_2DSp = Aerodynamic_forces(:,3);
        L = L_2DSp'*imp_points.W{z};

        D_ind2D_Sp = Aerodynamic_forces(:,1);
        D_i = D_ind2D_Sp'*imp_points.W{z};

        D_ind_2D=sum([reshape(D_ind2D_Sp, N_pannels_y(z), N_pannels_x(z))]');
        L_2D=sum([reshape(L_2DSp, N_pannels_y(z), N_pannels_x(z))]');

        induced_alpha=-asind(D_ind_2D./abs(L_2D));
    
        CD_2d = D_ind_2D./(norm(output_flux.U_inf{z})^2 *output_flux.rho * 0.5 * [imp_points.Choord{z}]' );
        CL_2d = L_2D./(norm(output_flux.U_inf{z})^2 *output_flux.rho * 0.5 * [imp_points.Choord{z}]' );

        CL = L/(norm(output_flux.U_inf{z})^2 * output_flux.rho * 0.5 * wing_data.Sup(1));
        CD = D_i/(norm(output_flux.U_inf{z})^2 * output_flux.rho * 0.5 * wing_data.Sup(1));
    
        sol.gamma{z} = gamma;

        sol.CL{z} = CL;
        sol.CD{z} = CD;

        sol.L{z} = L;
        sol.D_ind{z} = D_i;

        sol.ind_alpha{z}=induced_alpha;

        sol.D_ind2D{z}=D_ind_2D;  
        sol.L_2D{z}=L_2D;

        sol.Cd_2d{z}=CD_2d;
        sol.CL_2d{z}=CL_2d;
    
        i_p = N_pannels_body1;
    end
    
    clear L D_i CD CL gamma induced_alpha_Sp gamma_sp

    %% Noteworthy results
    
    L=0;
    D_i=0;
    CL=0;
    
    for z=1:config.NBodies
        L=L+sum(sol.L{z}(:,1));
        D_i=D_i+sum(sol.D_ind{z}(:));
        CL=CL+sol.CL{z};
    end

    config.CD=D_i/(U_inf_mod^2 * output_flux.rho * 0.5 * wing_data.Sup(1));
    config.CL=CL;
    
    %% PLOT

    if plot_flag==1

    fprintf('Total Lift = %.3f N\n',L);
    fprintf('Total induced drag = %.3f N\n',D_i);

        figure
        hold on
        for z = 1: config.NBodies
        
            X = plot_data.X{z};
            Y = plot_data.Y{z};
            Z = plot_data.Z{z};
        
            surf(X, Y, Z, 'FaceColor', 'cyan', 'EdgeColor', 'k');

            xlabel('X', Fontsize=12)
            ylabel('Y', Fontsize=12)
            zlabel('Z', Fontsize=12, Rotation=0)
            grid on;
            axis("image")
            view(3)
            title('Geometry', FontSize=18)

        end
        hold off
        
        %% plot circ su ala 
        
        figure 
        hold on
        for z = 1: config.NBodies
        
            gamma_sp=sol.gamma_sp{z};
        
            X = plot_data.X{z};
            Y = plot_data.Y{z};
            Z = plot_data.Z{z};
        
            [Nx_plot,Ny_plot] = size(gamma_sp);
        
            for j = 1:Ny_plot
                for i =  1:Nx_plot
                    gamma = gamma_sp(i,j); 

                    X_panel = [X(i, j), X(i+1, j), X(i+1, j+1), X(i, j+1)];
                    Y_panel = [Y(i, j), Y(i+1, j), Y(i+1, j+1), Y(i, j+1)];
                    Z_panel = [Z(i, j), Z(i+1, j), Z(i+1, j+1), Z(i, j+1)];

                    fill3(X_panel, Y_panel, Z_panel,gamma, 'FaceAlpha', 0.7, 'EdgeColor', 'k');
                end
            end
        
            xlabel('X', Fontsize=12)
            ylabel('Y', Fontsize=12)
            zlabel('Z', Fontsize=12, Rotation=0)
            grid on
            axis("image")
            view(3)   
            colormap jet   
            c=colorbar;     
            c.Label.String = '\fontsize{15} \Gamma';
            title('Circulation on each panel', FontSize=18)
        
        end
        
        %% plot lift/drag 

        figure
         
        for z = 1:config.NBodies
            hold on
            xx = linspace(-config.SemiSpan(z),config.SemiSpan(z),N_pannels_y(z));
            plot(xx, sol.L_2D{1,z},'o-','LineWidth',1.5);
            grid on 
            title('2D Lift distribution', FontSize=18);
            xlabel('Spanwise position [m]', FontSize=12);
            ylabel('L [N/m]', Rotation=0, FontSize=12)
        end

        legend({'Wing', 'Tail'}, 'Location', 'best')
        %%

        figure
        for z = 1:config.NBodies
            hold on
            xx2 = linspace(-config.SemiSpan(z),config.SemiSpan(z),N_pannels_y(z));
            plot(xx2,sol.ind_alpha{1,z},'o-','LineWidth',1.5);
            grid on 
            title('Induced angles', FontSize=18);
            xlabel('Spanwise position [m]', FontSize=12);
            ylabel('\alpha_i [deg]', Rotation=0, FontSize=12);
        end
        legend({'Wing', 'Tail'}, 'Location', 'best')

        figure
        for z = 1:config.NBodies
            % Visualize induced drag per units of length
            hold on
            xx2 = linspace(-config.SemiSpan(z),config.SemiSpan(z),N_pannels_y(z));
            plot(xx2,sol.D_ind2D{1,z} ,'o-','LineWidth',1.5);
            grid on 
            title('Induced Drag', FontSize=18);
            xlabel('Spanwise position [m]', 'FontSize',12);
            ylabel('D_i [N/m]', 'Rotation', 0, FontSize=12);
        end
        legend({'Wing', 'Tail'}, 'Location', 'best')

        clear xx xx2
    end

    Cl_zita(f, 1)=config.CL;
    Cd_zita(f, 1)=config.CD;

end