
function[imp_points,plot]=create_Geometry...
    (structure_data,wing_data,imp_points,plot,z)

if wing_data.d(z)~=0

    if wing_data.cmw(z)==wing_data.cr(z)

        ba = wing_data.b(z); 
        Croot = wing_data.cr(z); 
        Ctip = wing_data.ct(z); 
        d = wing_data.d(z);
        Sweep1 = wing_data.sweep1(z);
        Sweep2 = wing_data.sweep2(z);
        Dihedral = wing_data.diedro(z);
        
        N_pannels_x = structure_data.choord_panels(z);
        N_pannels_y = structure_data.semiSpan_panels(z);
        Pos_x = structure_data.LEposition_x(z);
        Pos_z = structure_data.LEposition_z(z);
        RotationA_y = structure_data.RotationA_y(z);
        
        l = (ba/2) - d;
 
        Sweepangle1 = deg2rad(Sweep1+90); 
        Sweepangle2 = deg2rad(Sweep2+90);
        Dihedral_angle = deg2rad(Dihedral); 
        Rotation_alfa = deg2rad(RotationA_y);     
        
        Cos_alfa = cos(Rotation_alfa);
        Sin_alfa = sin(Rotation_alfa);
        
        R_alfa = [Cos_alfa 0 Sin_alfa; 0 1 0; -Sin_alfa 0 Cos_alfa]; 
          
        Cos_Dihedral = cos(Dihedral_angle);
        Sin_Dihedral = sin(Dihedral_angle);
        Rest = rem(N_pannels_y, 2);  
        
        R_dihedral = [1 0 0; 0 Cos_Dihedral -Sin_Dihedral; 0 Sin_Dihedral Cos_Dihedral];
         
        if Rest == 0
            ny_rect_e = round(d/(ba/2)*N_pannels_y)+1; 
            ny_trap_e = round(l/(ba/2)*N_pannels_y)+1; 
        else
            ny_rect_e = floor(d/(ba/2)*N_pannels_y)+1;
            ny_trap_e = ceil(l/(ba/2)*N_pannels_y)+1;
        end
        
        N_pannels_x = N_pannels_x +1; 

        %% rectangular

        y1 = linspace(0,d,ny_rect_e);
        y_punti_centrali_1 = (y1(1:end-1) + y1(2:end)) / 2; 
        Ps1 = -d*cot(Sweepangle1); 
        
        if abs(Ps1)<1e-8
            Ps1=0; 
        end
        
        A1 = [0 , y1(1)];
        B1 = [Croot, y1(1)];
        D1 = [Ps1, y1(end)];
        C1 = [Ps1 + Croot, y1(end)];
        
        x_Le_1 = linspace(A1(1),D1(1),ny_rect_e);
        x_Ri_1 = linspace(B1(1),C1(1),ny_rect_e);
        
        %% quarter point
        
        X1_point_1 = zeros(N_pannels_x, ny_rect_e);  
        delta1_x_1 = zeros(N_pannels_x-1,ny_rect_e); 
        X_point_lat_1 = zeros(N_pannels_x-1,ny_rect_e);
        Y_point_1 = repmat(y1,N_pannels_x, 1);
        
        for j = 1:ny_rect_e
            X1_point_1(:,j) = linspace(x_Le_1(j), x_Ri_1(j), N_pannels_x); 
            
            for i = 1:N_pannels_x-1
                delta1_x_1(i,j) = X1_point_1(i+1,j) - X1_point_1(i,j);
                X_point_lat_1(i,j) = X1_point_1(i,j) + delta1_x_1(i,j)/4; 
            end
                  
        end
        
        Y_point_lat_1 = Y_point_1(1:end-1,:); 
        
        %% three quarter point 
        
        X1_point_2 = zeros(N_pannels_x,ny_rect_e-1);
        delta1_x_2 = zeros(N_pannels_x-2,ny_rect_e-1);
        X_point_cen_1 = zeros(N_pannels_x-2,ny_rect_e-1);
        Y_point_cen_1 = repmat(y_punti_centrali_1,N_pannels_x-1,1);
        
        for j = 1:ny_rect_e-1
            Y_point_cen_1(:,j) = y_punti_centrali_1(j) * ones(N_pannels_x-1,1); % xx
            for i = 1:N_pannels_x
                X1_point_2(i,j) = (X1_point_1(i,j+1) + X1_point_1(i,j))/2;
            end
        end
        
        for j = 1:ny_rect_e-1
            for i = 2:N_pannels_x
                delta1_x_2(i-1,j) = X1_point_2(i,j) - X1_point_2(i-1,j);
            end
        end
        
        for j = 1:ny_rect_e-1
            for i = 1:N_pannels_x-1
            X_point_cen_1(i,j)= X1_point_2(i,j) + 3*delta1_x_2(i,j) /4;
            end
        end
        
        %% trapezoidal
 
        y2 = linspace(d,ba/2,ny_trap_e);
        y_punti_centrali_2 = (y2(1:end-1) + y2(2:end)) / 2; 
        Ps2 = -l*cot(Sweepangle2);
        A2 = [D1(1),y1(end)];
        B2 = [C1(1),y1(end)];
        D2 = [D1(1)+Ps2,y2(end)];
        C2 = [D1(1)+Ps2+Ctip,y2(end)];
        
        x_left_2 = linspace(A2(1),D2(1),ny_trap_e);
        x_right_2 = linspace(B2(1),C2(1),ny_trap_e);
        
        %%  quarter point 2
        
        X2_point_1 = zeros(N_pannels_x, ny_trap_e);
        delta2_x_1 = zeros(N_pannels_x-1,ny_trap_e);
        X_point_lat_2 = zeros(N_pannels_x-1,ny_trap_e);
        Y_point_2 = repmat(y2,N_pannels_x, 1);
        
        for j = 1:ny_trap_e
            X2_point_1(:,j) = linspace(x_left_2(j), x_right_2(j), N_pannels_x);
            
            for i = 1:N_pannels_x-1
                delta2_x_1(i,j) = X2_point_1(i+1,j) - X2_point_1(i,j);
                X_point_lat_2(i,j) = X2_point_1(i,j) + delta2_x_1(i,j)/4;
            end
                  
        end
        
        Y_punti_laterali_2 = Y_point_2(1:end-1,:);
    
        %% three quarter point 2
        
        X2_point_2 = zeros(N_pannels_x,ny_trap_e-1);
        delta2_x_2 = zeros(N_pannels_x-2,ny_trap_e-1);
        X_point_cen_2 = zeros(N_pannels_x-2,ny_trap_e-1);
        Y_point_cen_2 = repmat(y_punti_centrali_2,N_pannels_x-1,1);
        
        for j = 1:ny_trap_e-1
            for i = 1:N_pannels_x
                X2_point_2(i,j) = (X2_point_1(i,j+1) + X2_point_1(i,j))/2;
            end
        end
        
        for j = 1:ny_trap_e-1
            for i = 2:N_pannels_x
                delta2_x_2(i-1,j) = X2_point_2(i,j) - X2_point_2(i-1,j);
            end
        end
        
        for j = 1:ny_trap_e-1
            for i = 1:N_pannels_x-1
            X_point_cen_2(i,j)= X2_point_2(i,j) + 3*delta2_x_2(i,j) /4;
            end
        end
        
        %% Rotation on Y 

        Cos_alfa = R_alfa(1,1);
        Sin_alfa = R_alfa(3,1);
        
        Z_punti_1 = zeros(size(Y_point_1));
        X1_point_1 = Cos_alfa * X1_point_1 - Sin_alfa * Z_punti_1;
        Z_punti_1 = Sin_alfa * X1_point_1 + Cos_alfa * Z_punti_1; 
        
        Z_punti_laterali_1 = zeros(size(Y_point_lat_1));
        X_point_lat_1 = Cos_alfa * X_point_lat_1 - Sin_alfa * Z_punti_laterali_1;
        Z_punti_laterali_1 = Sin_alfa * X_point_lat_1 + Cos_alfa * Z_punti_laterali_1; 
        
        Z_punti_centrali_1 = zeros(size(Y_point_cen_1));
        X_point_cen_1 = Cos_alfa * X_point_cen_1 - Sin_alfa * Z_punti_centrali_1;
        Z_punti_centrali_1 = Sin_alfa * X_point_cen_1 + Cos_alfa * Z_punti_centrali_1; 
        
        Z_punti_2 = zeros(size(Y_point_2));
        X2_point_1 = Cos_alfa * X2_point_1 - Sin_alfa * Z_punti_2;
        Z_punti_2 = Sin_alfa * X2_point_1 + Cos_alfa * Z_punti_2; 
        
        Z_punti_laterali_2 = zeros(size(Y_punti_laterali_2));
        X_point_lat_2 = Cos_alfa * X_point_lat_2 - Sin_alfa * Z_punti_laterali_2;
        Z_punti_laterali_2 = Sin_alfa * X_point_lat_2 + Cos_alfa * Z_punti_laterali_2; 
        
        Z_punti_centrali_2 = zeros(size(Y_point_cen_2));
        X_point_cen_2 = Cos_alfa * X_point_cen_2 - Sin_alfa * Z_punti_centrali_2;
        Z_punti_centrali_2 = Sin_alfa * X_point_cen_2 + Cos_alfa * Z_punti_centrali_2;
        
        %% Dihedral angle 
        
        R_a_d = R_alfa * R_dihedral * R_alfa';
        
        X1_point_1 = R_a_d(1,1) * X1_point_1 + R_a_d(1,2) * Y_point_1 + R_a_d(1,3)*Z_punti_1;
        Y_punti_1_rot = R_a_d(2,1) * X1_point_1 + R_a_d(2,2) * Y_point_1 + R_a_d(2,3)*Z_punti_1;
        Z_punti_1_rot = R_a_d(3,1) * X1_point_1 + R_a_d(3,2) * Y_point_1 + R_a_d(3,3)*Z_punti_1; 
        
        X_point_lat_1 = R_a_d(1,1) * X_point_lat_1 + R_a_d(1,2) * Y_point_lat_1 + R_a_d(1,3)*Z_punti_laterali_1;
        Y_punti_laterali_1_rot = R_a_d(2,1) * X_point_lat_1 + R_a_d(2,2) * Y_point_lat_1 + R_a_d(2,3)*Z_punti_laterali_1;
        Z_punti_laterali_1_rot = R_a_d(3,1) * X_point_lat_1 + R_a_d(3,2) * Y_point_lat_1 + R_a_d(3,3)*Z_punti_laterali_1; 
        
        X_point_cen_1 = R_a_d(1,1) * X_point_cen_1 + R_a_d(1,2) * Y_point_cen_1 + R_a_d(1,3)*Z_punti_centrali_1;
        Y_punti_centrali_1_rot = R_a_d(2,1) * X_point_cen_1 + R_a_d(2,2) * Y_point_cen_1 + R_a_d(2,3)*Z_punti_centrali_1;
        Z_punti_centrali_1_rot = R_a_d(3,1) * X_point_cen_1 + R_a_d(3,2) * Y_point_cen_1 + R_a_d(3,3)*Z_punti_centrali_1; 
        
        X2_point_1 = R_a_d(1,1) * X2_point_1 + R_a_d(1,2) * Y_point_2 + R_a_d(1,3)*Z_punti_2;
        Y_punti_2_rot =  R_a_d(2,1) * X2_point_1 + R_a_d(2,2) * Y_point_2 + R_a_d(2,3)*Z_punti_2;
        Z_punti_2_rot = R_a_d(3,1) * X2_point_1 + R_a_d(3,2) * Y_point_2 + R_a_d(3,3)*Z_punti_2;
        
        X_point_lat_2 = R_a_d(1,1) * X_point_lat_2 + R_a_d(1,2) * Y_punti_laterali_2 + R_a_d(1,3)*Z_punti_laterali_2;
        Y_punti_laterali_2_rot = R_a_d(2,1) * X_point_lat_2 + R_a_d(2,2) * Y_punti_laterali_2 + R_a_d(2,3)*Z_punti_laterali_2;
        Z_punti_laterali_2_rot = R_a_d(3,1) * X_point_lat_2 + R_a_d(3,2) * Y_punti_laterali_2 + R_a_d(3,3)*Z_punti_laterali_2; 
        
        X_point_cen_2 = R_a_d(1,1) * X_point_cen_2 + R_a_d(1,2) * Y_point_cen_2 + R_a_d(1,3)*Z_punti_centrali_2;
        Y_punti_centrali_2_rot = R_a_d(2,1) * X_point_cen_2 + R_a_d(2,2) * Y_point_cen_2 + R_a_d(2,3)*Z_punti_centrali_2;
        Z_punti_centrali_2_rot = R_a_d(3,1) * X_point_cen_2 + R_a_d(3,2) * Y_point_cen_2 + R_a_d(3,3)*Z_punti_centrali_2; 
        
        %% 
        
        Y_punti_1_rot_sx = -(Y_punti_1_rot);
        Z_punti_1_rot_sx = (Z_punti_1_rot);
        
        Y_punti_laterali_1_rot_sx = -(Y_punti_laterali_1_rot);
        Z_punti_laterali_1_rot_sx = (Z_punti_laterali_1_rot) ;
        
        Y_punti_centrali_1_rot_sx = -(Y_punti_centrali_1_rot);
        Z_punti_centrali_1_rot_sx = (Z_punti_centrali_1_rot);
        
        Y_punti_2_rot_sx = -(Y_punti_2_rot);
        Z_punti_2_rot_sx = (Z_punti_2_rot);
        
        Y_punti_laterali_2_rot_sx = -(Y_punti_laterali_2_rot);
        Z_punti_laterali_2_rot_sx = (Z_punti_laterali_2_rot);
        
        Y_punti_centrali_2_rot_sx = -(Y_punti_centrali_2_rot);
        Z_punti_centrali_2_rot_sx = (Z_punti_centrali_2_rot);
        
        %%
        
        X = Pos_x + [fliplr(X2_point_1(:,2:end)),fliplr(X1_point_1),X1_point_1(:,2:end),X2_point_1(:,2:end)];
        Y = [fliplr(Y_punti_2_rot(:,2:end)),fliplr(Y_punti_1_rot),Y_punti_1_rot_sx(:,2:end),Y_punti_2_rot_sx(:,2:end)];
        Z = Pos_z + [fliplr(Z_punti_2_rot(:,2:end)),fliplr(Z_punti_1_rot),Z_punti_1_rot_sx(:,2:end),Z_punti_2_rot_sx(:,2:end)];
    
        X_punti_laterali_1_sx = fliplr(X_point_lat_1);
        X_punti_laterali_2_sx = fliplr(X_point_lat_2);
        X_punti_centrali_1_sx = fliplr(X_point_cen_1);
        X_punti_centrali_2_sx = fliplr(X_point_cen_2);
        
        X_l_t = Pos_x + [X_punti_laterali_2_sx(:,1:end-1),X_punti_laterali_1_sx(:,1:end-1),X_point_lat_1,X_point_lat_2(:,2:end)];
        Y_l_t = [fliplr(Y_punti_laterali_2_rot_sx(:,2:end)),fliplr(Y_punti_laterali_1_rot_sx(:,2:end)),Y_punti_laterali_1_rot,Y_punti_laterali_2_rot(:,2:end)];
        Z_l_t = Pos_z + [fliplr(Z_punti_laterali_2_rot_sx(:,2:end)),fliplr(Z_punti_laterali_1_rot_sx(:,2:end)),Z_punti_laterali_1_rot,Z_punti_laterali_2_rot(:,2:end)];
        
        X_c_t = Pos_x + [X_punti_centrali_2_sx,X_punti_centrali_1_sx,X_point_cen_1,X_point_cen_2];
        Y_c_t= [fliplr(Y_punti_centrali_2_rot_sx),fliplr(Y_punti_centrali_1_rot_sx),Y_punti_centrali_1_rot,Y_punti_centrali_2_rot];
        Z_c_t = Pos_z + [fliplr(Z_punti_centrali_2_rot_sx),fliplr(Z_punti_centrali_1_rot_sx),Z_punti_centrali_1_rot,Z_punti_centrali_2_rot];
        
        plot.X{z} = X;
        plot.Y{z} = Y;
        plot.Z{z} = Z;
        
        plot.Ex{z} = X_l_t; 
        plot.Ey{z} = Y_l_t;
        plot.Ez{z} = Z_l_t;
        
        plot.Cx{z} = X_c_t;
        plot.Cy{z} = Y_c_t;
        plot.Cz{z} = Z_c_t;
        
        %%
        
        [nxt,nyt] = size(X_c_t); 
        Np = nyt*nxt;
        E1 = zeros(Np,3); 
        E2 = zeros(Np,3); 
        Center = zeros(Np,3);
        
        normdx = R_alfa * [0 ; 0 ; 1];
        normdx = R_a_d * normdx;
        normsx = normdx.*[1 ; -1 ; 1];
        normal_sx = [ones(nyt/2,1)*normsx(1),ones(nyt/2,1)*normsx(2), ones(nyt/2,1)*normsx(3)];
        normal_dx = [ones(nyt/2,1)*normdx(1),ones(nyt/2,1)*normdx(2), ones(nyt/2,1)*normdx(3)];
        Normal = zeros(Np,3);
    
        tangentdx=R_alfa * [1 ; 0 ; 0];
        tangentdx=R_a_d * tangentdx;
        tangentsx=tangentdx.*[1 ; -1 ; 1];
        tangent_sx=[ones(nyt/2,1)*tangentsx(1),ones(nyt/2,1)*tangentsx(2), ones(nyt/2,1)*tangentsx(3)];
        tangent_dx=[ones(nyt/2,1)*tangentdx(1),ones(nyt/2,1)*tangentdx(2), ones(nyt/2,1)*tangentdx(3)];
        Tangent=zeros(Np, 3);
    
        for i = 1: nxt
            ind1 = nyt*(i-1)+1;
            ind2 = i*nyt;
            E1(ind1:ind2,:) = [X_l_t(i,1:end-1)',Y_l_t(i,1:end-1)',Z_l_t(i,1:end-1)'];
            E2(ind1:ind2,:) = [X_l_t(i,2:end)',Y_l_t(i,2:end)',Z_l_t(i,2:end)'];
            Center(ind1:ind2,:) = [X_c_t(i,:)',Y_c_t(i,:)',Z_c_t(i,:)'];
            Normal(ind1:ind2,:) = [normal_sx;normal_dx];
            W = abs(E1(:,2)-E2(:,2)); 
            Tangent(ind1:ind2, :)=[tangent_sx; tangent_dx]; 
        end
    
        choor_trap1 = Croot*ones(ny_rect_e-1,1); % corde panelli tratto rettangolare 
       
        choord_trap2 = zeros(ny_trap_e-1,1);
        delta_c2 = (Croot - Ctip)/((ny_trap_e-1));
        choord_trap2(1,1) = Croot - (delta_c2/2);
        for i = 2: ny_trap_e-1
            choord_trap2(i,1) = choord_trap2(i-1,1) - delta_c2;
        end
    
        choord_sim_section = [choor_trap1;choord_trap2];
        Choord = [flipud(choord_sim_section);choord_sim_section];
        
        imp_points.E1{z} = E1;
        imp_points.E2{z} = E2;
        imp_points.Center{z} = Center;

        imp_points.W{z} = W;
        imp_points.Choord{z} = Choord;

        imp_points.Normal{z} = Normal;
        imp_points.Tangent{z}=Tangent;

    end

    %%
%
%
%
%
%
%
%

    if wing_data.cmw(z)~=wing_data.cr(z)

        ba = wing_data.b(z); 
        Croot = wing_data.cr(z); 
        Ctip = wing_data.ct(z); 
        d = wing_data.d(z); 
        CM = wing_data.cmw(z);
        Sweep1 = wing_data.sweep1(z);
        Sweep2 = wing_data.sweep2(z);
        Dihedral = wing_data.diedro(z);
        
        N_pannels_x = structure_data.choord_panels(z);
        N_pannels_y = structure_data.semiSpan_panels(z);
        Pos_x = structure_data.LEposition_x(z);
        Pos_z = structure_data.LEposition_z(z);
        RotationA_y = structure_data.RotationA_y(z);
        
        l = (ba/2) - d; 
        
        Sweepangle1 = deg2rad(Sweep1+90); 
        Sweepangle2 = deg2rad(Sweep2+90);
        Dihedral_angle = deg2rad(Dihedral); 
        Rotation_alfa = deg2rad(RotationA_y);    
        
        Cos_alfa = cos(Rotation_alfa);
        Sin_alfa = sin(Rotation_alfa);
        
        R_alfa = [Cos_alfa 0 Sin_alfa; 0 1 0; -Sin_alfa 0 Cos_alfa]; 
          
        Cos_Dihedral = cos(Dihedral_angle);
        Sin_Dihedral = sin(Dihedral_angle);
        Rest= rem(N_pannels_y, 2); % trova il resto della divisione 
        
        R_dihedral = [1 0 0; 0 Cos_Dihedral -Sin_Dihedral; 0 Sin_Dihedral Cos_Dihedral];
         
        if Rest == 0
            ny_rect_e = round(d/(ba/2)*N_pannels_y)+1; % numero estremi pannello lungo y rettangolo 
            ny_trap_e = round(l/(ba/2)*N_pannels_y)+1; % numero estremi pannello lungo y trapezio
        else
            ny_rect_e = floor(d/(ba/2)*N_pannels_y)+1; % stessa cosa di sopra ma approssima in altro modo 
            ny_trap_e = ceil(l/(ba/2)*N_pannels_y)+1;
        end
        
        N_pannels_x = N_pannels_x +1; %numero di discretizzazioni
        %% trapezoidal 1
        
        % Vertici del trapezio 
        y1 = linspace(0,d,ny_rect_e);
        y_punti_centrali_1 = (y1(1:end-1) + y1(2:end)) / 2; 
        Ps1 = -d*cot(Sweepangle1); 
        
        if abs(Ps1)<1e-8
            Ps1=0; 
        end
        
        A1 = [0 , y1(1)];
        B1 = [Croot, y1(1)];
        D1 = [Ps1, y1(end)];
        C1 = [Ps1 + CM, y1(end)];
    
        x_Le_1 = linspace(A1(1),D1(1),ny_rect_e);
        x_Ri_1 = linspace(B1(1),C1(1),ny_rect_e);
        
        % quarter point
        
        % calcolo estremi dei panelli e quarto di corda per ogni panello 
        
        X1_point_1 = zeros(N_pannels_x, ny_rect_e);  % punti estremi panelli lungo asse x 
        delta1_x_1 = zeros(N_pannels_x-1,ny_rect_e); % ampiezza dei panelli lungo x
        X_point_lat_1 = zeros(N_pannels_x-1,ny_rect_e);
        Y_point_1 = repmat(y1,N_pannels_x, 1);
        
        for j = 1:ny_rect_e
            X1_point_1(:,j) = linspace(x_Le_1(j), x_Ri_1(j), N_pannels_x); 
            
            for i = 1:N_pannels_x-1
                delta1_x_1(i,j) = X1_point_1(i+1,j) - X1_point_1(i,j);
                X_point_lat_1(i,j) = X1_point_1(i,j) + delta1_x_1(i,j)/4;  % altezza a 1/4 ddella corda per ogni panello 
            end
                  
        end
        
        Y_point_lat_1 = Y_point_1(1:end-1,:); % estremità di ogni panello  
        
        % three quarter point 
        
        X1_point_2 = zeros(N_pannels_x,ny_rect_e-1);
        delta1_x_2 = zeros(N_pannels_x-2,ny_rect_e-1);
        X_point_cen_1 = zeros(N_pannels_x-2,ny_rect_e-1);
        Y_point_cen_1 = repmat(y_punti_centrali_1,N_pannels_x-1,1);
        
        for j = 1:ny_rect_e-1
            for i = 1:N_pannels_x
                X1_point_2(i,j) = (X1_point_1(i,j+1) + X1_point_1(i,j))/2;
            end
        end
        
        for j = 1:ny_rect_e-1
            for i = 2:N_pannels_x
                delta1_x_2(i-1,j) = X1_point_2(i,j) - X1_point_2(i-1,j);
            end
        end
        
        for j = 1:ny_rect_e-1
            for i = 1:N_pannels_x-1
            X_point_cen_1(i,j)= X1_point_2(i,j) + 3*delta1_x_2(i,j) /4;
            end
        end
        
        %% trapezoidal 2
         
        y2 = linspace(d,ba/2,ny_trap_e);
        y_punti_centrali_2 = (y2(1:end-1) + y2(2:end)) / 2; 
        Ps2 = -l*cot(Sweepangle2);
        A2 = [D1(1),y1(end)];
        B2 = [C1(1),y1(end)];
        D2 = [D1(1)+Ps2,y2(end)];
        C2 = [D1(1)+Ps2+Ctip,y2(end)];
        
        x_left_2 = linspace(A2(1),D2(1),ny_trap_e);
        x_right_2 = linspace(B2(1),C2(1),ny_trap_e);
        
        %  quarter point 2
        
        % sto calcolando i punti estremi di ogni panello 
        
        X2_point_1 = zeros(N_pannels_x, ny_trap_e);
        delta2_x_1 = zeros(N_pannels_x-1,ny_trap_e);
        X_point_lat_2 = zeros(N_pannels_x-1,ny_trap_e);
        Y_point_2 = repmat(y2,N_pannels_x, 1);
        
        for j = 1:ny_trap_e
            X2_point_1(:,j) = linspace(x_left_2(j), x_right_2(j), N_pannels_x);
            
            for i = 1:N_pannels_x-1
                delta2_x_1(i,j) = X2_point_1(i+1,j) - X2_point_1(i,j);
                X_point_lat_2(i,j) = X2_point_1(i,j) + delta2_x_1(i,j)/4;
            end
                  
        end
        
        Y_punti_laterali_2 = Y_point_2(1:end-1,:);
    
        % three quarter point 2
        
        X2_point_2 = zeros(N_pannels_x,ny_trap_e-1);
        delta2_x_2 = zeros(N_pannels_x-2,ny_trap_e-1);
        X_point_cen_2 = zeros(N_pannels_x-2,ny_trap_e-1);
        Y_point_cen_2 = repmat(y_punti_centrali_2,N_pannels_x-1,1);
        
        for j = 1:ny_trap_e-1
            Y_point_cen_2(:,j) = y_punti_centrali_2(j) * ones(N_pannels_x-1,1); % xx
            for i = 1:N_pannels_x
                X2_point_2(i,j) = (X2_point_1(i,j+1) + X2_point_1(i,j))/2;
            end
        end
        
        for j = 1:ny_trap_e-1
            for i = 2:N_pannels_x
                delta2_x_2(i-1,j) = X2_point_2(i,j) - X2_point_2(i-1,j);
            end
        end
        
        for j = 1:ny_trap_e-1
            for i = 1:N_pannels_x-1
            X_point_cen_2(i,j)= X2_point_2(i,j) + 3*delta2_x_2(i,j) /4;
            end
        end
        
        % fino a qui controllato riga per riga 
        
        %% rotation on Y 

        Cos_alfa = R_alfa(1,1);
        Sin_alfa = R_alfa(3,1);
        
        Z_punti_1 = zeros(size(Y_point_1));
        X1_point_1 = Cos_alfa * X1_point_1 - Sin_alfa * Z_punti_1;
        Z_punti_1 = Sin_alfa * X1_point_1 + Cos_alfa * Z_punti_1; 
        
        Z_punti_laterali_1 = zeros(size(Y_point_lat_1));
        X_point_lat_1 = Cos_alfa * X_point_lat_1 - Sin_alfa * Z_punti_laterali_1;
        Z_punti_laterali_1 = Sin_alfa * X_point_lat_1 + Cos_alfa * Z_punti_laterali_1; 
        
        Z_punti_centrali_1 = zeros(size(Y_point_cen_1));
        X_point_cen_1 = Cos_alfa * X_point_cen_1 - Sin_alfa * Z_punti_centrali_1;
        Z_punti_centrali_1 = Sin_alfa * X_point_cen_1 + Cos_alfa * Z_punti_centrali_1; 
        
        Z_punti_2 = zeros(size(Y_point_2));
        X2_point_1 = Cos_alfa * X2_point_1 - Sin_alfa * Z_punti_2;
        Z_punti_2 = Sin_alfa * X2_point_1 + Cos_alfa * Z_punti_2; 
        
        Z_punti_laterali_2 = zeros(size(Y_punti_laterali_2));
        X_point_lat_2 = Cos_alfa * X_point_lat_2 - Sin_alfa * Z_punti_laterali_2;
        Z_punti_laterali_2 = Sin_alfa * X_point_lat_2 + Cos_alfa * Z_punti_laterali_2; 
        
        Z_punti_centrali_2 = zeros(size(Y_point_cen_2));
        X_point_cen_2 = Cos_alfa * X_point_cen_2 - Sin_alfa * Z_punti_centrali_2;
        Z_punti_centrali_2 = Sin_alfa * X_point_cen_2 + Cos_alfa * Z_punti_centrali_2;
        
        %% application of Dihedral angle
        
        R_a_d = R_alfa * R_dihedral * R_alfa';
        
        X1_point_1 = R_a_d(1,1) * X1_point_1 + R_a_d(1,2) * Y_point_1 + R_a_d(1,3)*Z_punti_1;
        Y_punti_1_rot = R_a_d(2,1) * X1_point_1 + R_a_d(2,2) * Y_point_1 + R_a_d(2,3)*Z_punti_1;
        Z_punti_1_rot = R_a_d(3,1) * X1_point_1 + R_a_d(3,2) * Y_point_1 + R_a_d(3,3)*Z_punti_1; 
        
        X_point_lat_1 = R_a_d(1,1) * X_point_lat_1 + R_a_d(1,2) * Y_point_lat_1 + R_a_d(1,3)*Z_punti_laterali_1;
        Y_punti_laterali_1_rot = R_a_d(2,1) * X_point_lat_1 + R_a_d(2,2) * Y_point_lat_1 + R_a_d(2,3)*Z_punti_laterali_1;
        Z_punti_laterali_1_rot = R_a_d(3,1) * X_point_lat_1 + R_a_d(3,2) * Y_point_lat_1 + R_a_d(3,3)*Z_punti_laterali_1; 
        
        X_point_cen_1 = R_a_d(1,1) * X_point_cen_1 + R_a_d(1,2) * Y_point_cen_1 + R_a_d(1,3)*Z_punti_centrali_1;
        Y_punti_centrali_1_rot = R_a_d(2,1) * X_point_cen_1 + R_a_d(2,2) * Y_point_cen_1 + R_a_d(2,3)*Z_punti_centrali_1;
        Z_punti_centrali_1_rot = R_a_d(3,1) * X_point_cen_1 + R_a_d(3,2) * Y_point_cen_1 + R_a_d(3,3)*Z_punti_centrali_1; 
        
        X2_point_1 = R_a_d(1,1) * X2_point_1 + R_a_d(1,2) * Y_point_2 + R_a_d(1,3)*Z_punti_2;
        Y_punti_2_rot =  R_a_d(2,1) * X2_point_1 + R_a_d(2,2) * Y_point_2 + R_a_d(2,3)*Z_punti_2;
        Z_punti_2_rot = R_a_d(3,1) * X2_point_1 + R_a_d(3,2) * Y_point_2 + R_a_d(3,3)*Z_punti_2;
        
        X_point_lat_2 = R_a_d(1,1) * X_point_lat_2 + R_a_d(1,2) * Y_punti_laterali_2 + R_a_d(1,3)*Z_punti_laterali_2;
        Y_punti_laterali_2_rot = R_a_d(2,1) * X_point_lat_2 + R_a_d(2,2) * Y_punti_laterali_2 + R_a_d(2,3)*Z_punti_laterali_2;
        Z_punti_laterali_2_rot = R_a_d(3,1) * X_point_lat_2 + R_a_d(3,2) * Y_punti_laterali_2 + R_a_d(3,3)*Z_punti_laterali_2; 
        
        X_point_cen_2 = R_a_d(1,1) * X_point_cen_2 + R_a_d(1,2) * Y_point_cen_2 + R_a_d(1,3)*Z_punti_centrali_2;
        Y_punti_centrali_2_rot = R_a_d(2,1) * X_point_cen_2 + R_a_d(2,2) * Y_point_cen_2 + R_a_d(2,3)*Z_punti_centrali_2;
        Z_punti_centrali_2_rot = R_a_d(3,1) * X_point_cen_2 + R_a_d(3,2) * Y_point_cen_2 + R_a_d(3,3)*Z_punti_centrali_2; 
        
        %%  left semi-wing
        
        Y_punti_1_rot_sx = -(Y_punti_1_rot);
        Z_punti_1_rot_sx = (Z_punti_1_rot);
        
        Y_punti_laterali_1_rot_sx = -(Y_punti_laterali_1_rot);
        Z_punti_laterali_1_rot_sx = (Z_punti_laterali_1_rot) ;
        
        Y_punti_centrali_1_rot_sx = -(Y_punti_centrali_1_rot);
        Z_punti_centrali_1_rot_sx = (Z_punti_centrali_1_rot);
        
        Y_punti_2_rot_sx = -(Y_punti_2_rot);
        Z_punti_2_rot_sx = (Z_punti_2_rot);
        
        Y_punti_laterali_2_rot_sx = -(Y_punti_laterali_2_rot);
        Z_punti_laterali_2_rot_sx = (Z_punti_laterali_2_rot);
        
        Y_punti_centrali_2_rot_sx = -(Y_punti_centrali_2_rot);
        Z_punti_centrali_2_rot_sx = (Z_punti_centrali_2_rot);
        
        %% 
        
        X = Pos_x + [fliplr(X2_point_1(:,2:end)),fliplr(X1_point_1),X1_point_1(:,2:end),X2_point_1(:,2:end)];
        Y = [fliplr(Y_punti_2_rot(:,2:end)),fliplr(Y_punti_1_rot),Y_punti_1_rot_sx(:,2:end),Y_punti_2_rot_sx(:,2:end)];
        Z = Pos_z + [fliplr(Z_punti_2_rot(:,2:end)),fliplr(Z_punti_1_rot),Z_punti_1_rot_sx(:,2:end),Z_punti_2_rot_sx(:,2:end)];
    
        X_punti_laterali_1_sx = fliplr(X_point_lat_1);
        X_punti_laterali_2_sx = fliplr(X_point_lat_2);
        X_punti_centrali_1_sx = fliplr(X_point_cen_1);
        X_punti_centrali_2_sx = fliplr(X_point_cen_2);
        
        X_l_t = Pos_x + [X_punti_laterali_2_sx(:,1:end-1),X_punti_laterali_1_sx(:,1:end-1),X_point_lat_1,X_point_lat_2(:,2:end)];
        Y_l_t = [fliplr(Y_punti_laterali_2_rot_sx(:,2:end)),fliplr(Y_punti_laterali_1_rot_sx(:,2:end)),Y_punti_laterali_1_rot,Y_punti_laterali_2_rot(:,2:end)];
        Z_l_t = Pos_z + [fliplr(Z_punti_laterali_2_rot_sx(:,2:end)),fliplr(Z_punti_laterali_1_rot_sx(:,2:end)),Z_punti_laterali_1_rot,Z_punti_laterali_2_rot(:,2:end)];
        
        X_c_t = Pos_x + [X_punti_centrali_2_sx,X_punti_centrali_1_sx,X_point_cen_1,X_point_cen_2];
        Y_c_t= [fliplr(Y_punti_centrali_2_rot_sx),fliplr(Y_punti_centrali_1_rot_sx),Y_punti_centrali_1_rot,Y_punti_centrali_2_rot];
        Z_c_t = Pos_z + [fliplr(Z_punti_centrali_2_rot_sx),fliplr(Z_punti_centrali_1_rot_sx),Z_punti_centrali_1_rot,Z_punti_centrali_2_rot];
        
        plot.X{z} = X;
        plot.Y{z} = Y;
        plot.Z{z} = Z;
        
        plot.Ex{z} = X_l_t; 
        plot.Ey{z} = Y_l_t;
        plot.Ez{z} = Z_l_t;
        
        plot.Cx{z} = X_c_t; 
        plot.Cy{z} = Y_c_t;
        plot.Cz{z} = Z_c_t;
        
        %%
        
        [nxt,nyt] = size(X_c_t); 
        Np = nyt*nxt;
        E1 = zeros(Np,3); 
        E2 = zeros(Np,3); 
        Center = zeros(Np,3);
        
        normdx = R_alfa * [0 ; 0 ; 1];
        normdx = R_a_d * normdx;
        normsx = normdx.*[1 ; -1 ; 1];
        normal_sx = [ones(nyt/2,1)*normsx(1),ones(nyt/2,1)*normsx(2), ones(nyt/2,1)*normsx(3)];
        normal_dx = [ones(nyt/2,1)*normdx(1),ones(nyt/2,1)*normdx(2), ones(nyt/2,1)*normdx(3)];
        Normal = zeros(Np,3);
    
        tangentdx=R_alfa * [1 ; 0 ; 0];
        tangentdx=R_a_d * tangentdx;
        tangentsx=tangentdx.*[1 ; -1 ; 1];
        tangent_sx=[ones(nyt/2,1)*tangentsx(1),ones(nyt/2,1)*tangentsx(2), ones(nyt/2,1)*tangentsx(3)];
        tangent_dx=[ones(nyt/2,1)*tangentdx(1),ones(nyt/2,1)*tangentdx(2), ones(nyt/2,1)*tangentdx(3)];
        Tangent=zeros(Np, 3);
    
        for i = 1: nxt
            ind1 = nyt*(i-1)+1;
            ind2 = i*nyt;
            E1(ind1:ind2,:) = [X_l_t(i,1:end-1)',Y_l_t(i,1:end-1)',Z_l_t(i,1:end-1)'];
            E2(ind1:ind2,:) = [X_l_t(i,2:end)',Y_l_t(i,2:end)',Z_l_t(i,2:end)'];
            Center(ind1:ind2,:) = [X_c_t(i,:)',Y_c_t(i,:)',Z_c_t(i,:)'];
            Normal(ind1:ind2,:) = [normal_sx;normal_dx];
            W = abs(E1(:,2)-E2(:,2)); 
            Tangent(ind1:ind2, :)=[tangent_sx; tangent_dx]; 
        end
    
        choor_trap1 = zeros(ny_rect_e-1,1); 
        delta_c1= (Croot - CM)/((ny_rect_e-1));
        choor_trap1(1,1) = Croot - delta_c1;

        choord_trap2 = zeros(ny_trap_e-1,1);
        delta_c2 = (CM - Ctip)/((ny_trap_e-1));
        choord_trap2(1,1) = CM - (delta_c2/2);

        for i = 2:ny_rect_e-1
            choor_trap1(i,1) = choor_trap1(i-1,1) - delta_c1;
        end

        for i = 2: ny_trap_e-1
            choord_trap2(i,1) = choord_trap2(i-1,1) - delta_c2;
        end
    
        choord_sim_section = [choor_trap1;choord_trap2];
        Choord = [flipud(choord_sim_section);choord_sim_section];
        
        imp_points.E1{z} = E1;
        imp_points.E2{z} = E2;
        imp_points.Center{z} = Center;

        imp_points.Choord{z} = Choord;
        imp_points.W{z} = W;

        imp_points.Normal{z} = Normal;
        imp_points.Tangent{z}=Tangent;

    end

end

%%
%
%
%
%
%
if wing_data.d(z)==0

    Croot = wing_data.cr(z);  
    Ctip = wing_data.ct(z);  
    b=wing_data.b(z); 
    Sweep1 = wing_data.sweep1(z);
    Dihedral = wing_data.diedro(z);
    
    N_pannels_x = structure_data.choord_panels(z);
    N_pannels_y = structure_data.semiSpan_panels(z);
    Pos_x = structure_data.LEposition_x(z);
    Pos_z = structure_data.LEposition_z(z);
    RotationA_y = structure_data.RotationA_y(z);
   
    semib=b/2;
    Sweepangle1 = deg2rad(Sweep1+90); 
    Dihedral_angle = deg2rad(Dihedral); 
    Rotation_alfa = deg2rad(RotationA_y);   
    
    Cos_alfa = cos(Rotation_alfa);
    Sin_alfa = sin(Rotation_alfa);
    
    R_alfa = [Cos_alfa 0 Sin_alfa; 0 1 0; -Sin_alfa 0 Cos_alfa]; 
    
    
    Cos_Dihedral = cos(Dihedral_angle);
    Sin_Dihedral = sin(Dihedral_angle);
    
    R_dihedral = [1 0 0; 0 Cos_Dihedral -Sin_Dihedral; 0 Sin_Dihedral Cos_Dihedral];
    
    N_pannels_y = N_pannels_y +1;
    N_pannels_x = N_pannels_x +1;  
    %% semi delta wing
    
    y1 = linspace(0,semib,N_pannels_y);
    y_punti_centrali_1 = (y1(1:end-1) + y1(2:end)) / 2; 
    Ps1 = -semib*cot(Sweepangle1); 
    
    if abs(Ps1)<1e-8
        Ps1=0; 
    end
    
    A1 = [0 , y1(1)];
    B1 = [Croot, y1(1)];
    D1 = [Ps1, y1(end)];
    C1 = [Ps1 + Ctip, y1(end)];
 
    x_Le_1 = linspace(A1(1),D1(1),N_pannels_y);
    x_Ri_1 = linspace(B1(1),C1(1),N_pannels_y);
    
    %% quarter point 
    
    X1_point_1 = zeros(N_pannels_x, N_pannels_y);  
    delta1_x_1 = zeros(N_pannels_x-1,N_pannels_y); 
    X_point_lat_1 = zeros(N_pannels_x-1,N_pannels_y);
    Y_point_1 = repmat(y1,N_pannels_x, 1); 
    
    for j = 1:N_pannels_y
        X1_point_1(:,j) = linspace(x_Le_1(j), x_Ri_1(j), N_pannels_x); 
        
        for i = 1:N_pannels_x-1
            delta1_x_1(i,j) = X1_point_1(i+1,j) - X1_point_1(i,j);
            X_point_lat_1(i,j) = X1_point_1(i,j) + delta1_x_1(i,j)/4; 
        end
              
    end
    
    Y_point_lat_1 = Y_point_1(1:end-1,:);
    
    
    %% three quarter point 
    
    X1_point_2 = zeros(N_pannels_x,N_pannels_y-1);
    delta1_x_2 = zeros(N_pannels_x-2,N_pannels_y-1);
    X_point_cen_1 = zeros(N_pannels_x-2,N_pannels_y-1);
    Y_point_cen_1 = repmat(y_punti_centrali_1,N_pannels_x-1,1);
    
    for j = 1:N_pannels_y-1
        
        for i = 1:N_pannels_x
            X1_point_2(i,j) = (X1_point_1(i,j+1) + X1_point_1(i,j))/2;
        end
    end
    
    for j = 1:N_pannels_y-1
        for i = 2:N_pannels_x
            delta1_x_2(i-1,j) = X1_point_2(i,j) - X1_point_2(i-1,j);
        end
    end
    
    for j = 1:N_pannels_y-1
        for i = 1:N_pannels_x-1
        X_point_cen_1(i,j)= X1_point_2(i,j) + 3*delta1_x_2(i,j) /4;
        end
    end
    
    %% Rotation on Y
    Cos_alfa = R_alfa(1,1);
    Sin_alfa = R_alfa(3,1);
    
    Z_punti_1 = zeros(size(Y_point_1));
    X1_point_1 = Cos_alfa * X1_point_1 - Sin_alfa * Z_punti_1;
    Z_punti_1 = Sin_alfa * X1_point_1 + Cos_alfa * Z_punti_1; 
    
    Z_punti_laterali_1 = zeros(size(Y_point_lat_1));
    X_point_lat_1 = Cos_alfa * X_point_lat_1 - Sin_alfa * Z_punti_laterali_1;
    Z_punti_laterali_1 = Sin_alfa * X_point_lat_1 + Cos_alfa * Z_punti_laterali_1; 
    
    Z_punti_centrali_1 = zeros(size(Y_point_cen_1));
    X_point_cen_1 = Cos_alfa * X_point_cen_1 - Sin_alfa * Z_punti_centrali_1;
    Z_punti_centrali_1 = Sin_alfa * X_point_cen_1 + Cos_alfa * Z_punti_centrali_1;    
    
    %% Dihedral angle 

    R_a_d = R_alfa * R_dihedral * R_alfa';
    
    X1_point_1 = R_a_d(1,1) * X1_point_1 + R_a_d(1,2) * Y_point_1 + R_a_d(1,3)*Z_punti_1;
    Y_punti_1_rot = R_a_d(2,1) * X1_point_1 + R_a_d(2,2) * Y_point_1 + R_a_d(2,3)*Z_punti_1;
    Z_punti_1_rot = R_a_d(3,1) * X1_point_1 + R_a_d(3,2) * Y_point_1 + R_a_d(3,3)*Z_punti_1; 
    
    X_point_lat_1 = R_a_d(1,1) * X_point_lat_1 + R_a_d(1,2) * Y_point_lat_1 + R_a_d(1,3)*Z_punti_laterali_1;
    Y_punti_laterali_1_rot = R_a_d(2,1) * X_point_lat_1 + R_a_d(2,2) * Y_point_lat_1 + R_a_d(2,3)*Z_punti_laterali_1;
    Z_punti_laterali_1_rot = R_a_d(3,1) * X_point_lat_1 + R_a_d(3,2) * Y_point_lat_1 + R_a_d(3,3)*Z_punti_laterali_1; 
    
    X_point_cen_1 = R_a_d(1,1) * X_point_cen_1 + R_a_d(1,2) * Y_point_cen_1 + R_a_d(1,3)*Z_punti_centrali_1;
    Y_punti_centrali_1_rot = R_a_d(2,1) * X_point_cen_1 + R_a_d(2,2) * Y_point_cen_1 + R_a_d(2,3)*Z_punti_centrali_1;
    Z_punti_centrali_1_rot = R_a_d(3,1) * X_point_cen_1 + R_a_d(3,2) * Y_point_cen_1 + R_a_d(3,3)*Z_punti_centrali_1; 

    %%  
    
    Y_punti_1_rot_sx = -(Y_punti_1_rot);
    Z_punti_1_rot_sx = (Z_punti_1_rot);
    
    Y_punti_laterali_1_rot_sx = -(Y_punti_laterali_1_rot);
    Z_punti_laterali_1_rot_sx = (Z_punti_laterali_1_rot) ;
    
    Y_punti_centrali_1_rot_sx = -(Y_punti_centrali_1_rot);
    Z_punti_centrali_1_rot_sx = (Z_punti_centrali_1_rot);
    
    %% 
    
    X = Pos_x + [fliplr(X1_point_1),X1_point_1(:,2:end)];
    Y = [fliplr(Y_punti_1_rot),Y_punti_1_rot_sx(:,2:end)];
    Z = Pos_z + [fliplr(Z_punti_1_rot),Z_punti_1_rot_sx(:,2:end)];   
    
    X_punti_laterali_1_sx = fliplr(X_point_lat_1);
    X_punti_centrali_1_sx = fliplr(X_point_cen_1);
    
    X_l_t = Pos_x + [X_punti_laterali_1_sx(:,1:end-1),X_point_lat_1];
    Y_l_t = [fliplr(Y_punti_laterali_1_rot_sx(:,2:end)),Y_punti_laterali_1_rot];
    Z_l_t = Pos_z + [fliplr(Z_punti_laterali_1_rot_sx(:,2:end)),Z_punti_laterali_1_rot];
    
    X_c_t = Pos_x + [X_punti_centrali_1_sx,X_point_cen_1];
    Y_c_t= [fliplr(Y_punti_centrali_1_rot_sx),Y_punti_centrali_1_rot];
    Z_c_t = Pos_z + [fliplr(Z_punti_centrali_1_rot_sx),Z_punti_centrali_1_rot];
        
    plot.X{z} = X;
    plot.Y{z} = Y;
    plot.Z{z} = Z;
    
    plot.Ex{z} = X_l_t; 
    plot.Ey{z} = Y_l_t;
    plot.Ez{z} = Z_l_t;
    
    plot.Cx{z} = X_c_t; 
    plot.Cy{z} = Y_c_t;
    plot.Cz{z} = Z_c_t;
    
    
    %%
    
    [nxt,nyt] = size(X_c_t); 
    Np = nyt*nxt;
    E1 = zeros(Np,3); 
    E2 = zeros(Np,3); 
    Center = zeros(Np,3);
    
    normdx = R_alfa * [0 ; 0 ; 1];
    normdx = R_a_d * normdx;
    normsx = normdx.*[1 ; -1 ; 1];
    normal_sx = [ones(nyt/2,1)*normsx(1),ones(nyt/2,1)*normsx(2), ones(nyt/2,1)*normsx(3)];
    normal_dx = [ones(nyt/2,1)*normdx(1),ones(nyt/2,1)*normdx(2), ones(nyt/2,1)*normdx(3)];
    Normal = zeros(Np,3);  

    tangentdx=R_alfa * [1 ; 0 ; 0];
    tangentdx=R_a_d * tangentdx;
    tangentsx=tangentdx.*[1 ; -1 ; 1];
    tangent_sx=[ones(nyt/2,1)*tangentsx(1),ones(nyt/2,1)*tangentsx(2), ones(nyt/2,1)*tangentsx(3)];
    tangent_dx=[ones(nyt/2,1)*tangentdx(1),ones(nyt/2,1)*tangentdx(2), ones(nyt/2,1)*tangentdx(3)];
    Tangent=zeros(Np, 3);
    
    for i = 1: nxt
    
        ind1 = nyt*(i-1)+1;
        ind2 = i*nyt;
        E1(ind1:ind2,:) = [X_l_t(i,1:end-1)',Y_l_t(i,1:end-1)',Z_l_t(i,1:end-1)'];
        E2(ind1:ind2,:) = [X_l_t(i,2:end)',Y_l_t(i,2:end)',Z_l_t(i,2:end)'];
        Center(ind1:ind2,:) = [X_c_t(i,:)',Y_c_t(i,:)',Z_c_t(i,:)'];
        Normal(ind1:ind2,:) = [normal_sx;normal_dx];
        W = abs(E1(:,2)-E2(:,2)); 
        Tangent(ind1:ind2, :)=[tangent_sx; tangent_dx]; 
    end
    
    choor_trap1 = zeros(N_pannels_y-1,1);
    delta_c1 = (Croot - Ctip)/((N_pannels_y-1));
    choor_trap1(1,1) = Croot - (delta_c1/2);
    for i = 2: N_pannels_y-1
        choor_trap1(i,1) = choor_trap1(i-1,1) - delta_c1;
    end
    
    choord_sim_section = choor_trap1; 
    Choord = [flipud(choord_sim_section);choord_sim_section];
    
    imp_points.E1{z} = E1;
    imp_points.E2{z} = E2;
    imp_points.Center{z} = Center;

    imp_points.W{z} = W;
    imp_points.Choord{z} = Choord;

    imp_points.Normal{z} = Normal;  
    imp_points.Tangent{z}=Tangent;
end
end
