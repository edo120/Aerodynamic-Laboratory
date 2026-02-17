clc; clear; close all;

%% ----------------------------------------------------
% CASE DEFINITIONS
%% ----------------------------------------------------
cases = {
    'bl_E168_N10_a4.txt',  'E168', '10';
    'bl_E168_N3_a4.txt',   'E168', '3';
    'bl_E171_N10_a4.txt',  'E171', '10';
    'bl_E171_N3_a4.txt',   'E171', '3';
    'bl_0008_N10_a4.txt',  'NACA 0008', '10';
    'bl_0008_N3_a4.txt',   'NACA 0008', '3';
    'bl_0012_N10_a4.txt',  'NACA 0012', '10';
    'bl_0012_N3_a4.txt',   'NACA 0012', '3';
};

unique_airfoils = unique(cases(:,2), 'stable');
colors_master = lines(size(cases,1)); 

fprintf('\n--- Boundary Layer & Pressure Diagnostics ---\n')

for a = 1:length(unique_airfoils)
    current_airfoil = unique_airfoils{a};
    figure('Color', 'w', 'Units', 'pixels', 'Position', [100, 100, 1000, 850], 'Name', current_airfoil);
    
    % Subplot 1: Shape Factor (H)
    axH = subplot(2,1,1); hold on; grid on; box on;
    title(['Airfoil: ', current_airfoil, ' - Shape Factor (H)'], 'FontSize', 12);
    ylabel('H', 'FontWeight', 'bold');
    
    % Subplot 2: Pressure Coefficient (Cp)
    axP = subplot(2,1,2); hold on; grid on; box on;
    set(axP, 'YDir', 'reverse'); % Cp standard: negative up
    title('Pressure Coefficient (C_p) & Plateau Detection', 'FontSize', 12);
    xlabel('Chordwise Position (x/c)', 'FontWeight', 'bold');
    ylabel('C_p (Inverted)', 'FontWeight', 'bold');

    case_indices = find(strcmp(cases(:,2), current_airfoil));
    
    for idx_case = 1:length(case_indices)
        i = case_indices(idx_case);
        filename = cases{i,1};
        ncrit = cases{i,3};
        
        if ~exist(filename, 'file'), continue; end
        data = readmatrix(filename, 'FileType', 'text');
        
        % Data mapping
        x_raw  = data(:,2);
        y_raw  = data(:,3);
        Ue_raw = data(:,4); 
        Cf_raw = data(:,7);
        H_raw  = data(:,8);
        
        surfaces = {'upper', 'lower'};
        for s = 1:2
            idx = (s==1 & y_raw >= 0) | (s==2 & y_raw < 0);
            if ~any(idx), continue; end
            
            [x, sIdx] = sort(x_raw(idx));
            H = H_raw(idx); H = H(sIdx);
            Cf = Cf_raw(idx); Cf = Cf(sIdx);
            Ue = Ue_raw(idx); Ue = Ue(sIdx);
            Cp = 1 - Ue.^2; 
            
            lStyle = '-'; if s==2, lStyle = '--'; end
            surfLab = 'Up'; if s==2, surfLab = 'Lo'; end
            col = colors_master(i,:);
            fullLabel = sprintf('%s (N=%s) %s', current_airfoil, ncrit, surfLab);

            % --- Plot H ---
            plot(axH, x, H, lStyle, 'LineWidth', 1.5, 'Color', col, 'DisplayName', fullLabel);
            
            % --- Plot Cp ---
            plot(axP, x, Cp, lStyle, 'LineWidth', 1.5, 'Color', col, 'HandleVisibility', 'off');

            % --- LSB DETECTION ---
            sep_idx = find(Cf <= 1e-6, 1);
            if ~isempty(sep_idx)
                reattach_idx = find(H(sep_idx:end) < 1.8, 1) + sep_idx - 1;
                
                if ~isempty(reattach_idx)
                    x_sep = x(sep_idx);
                    x_rea = x(reattach_idx);
                    
                    % 1. Plot Markers on H Plot
                    plot(axH, x_sep, H(sep_idx), 'rs', 'MarkerFaceColor', 'r', 'HandleVisibility', 'off');
                    plot(axH, x_rea, H(reattach_idx), 'ko', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
                    
                    % 2. Identify and Plot Plateau on Cp Plot
                    % Plateau end: look for where dCp/dx starts to increase significantly
                    dCpdx = gradient(Cp, x);
                    plat_end_idx = find(x > x_sep & dCpdx > 0.5, 1); 
                    
                    if ~isempty(plat_end_idx)
                        x_plat_end = x(plat_end_idx);
                        % Using a thick line with no AlphaData to ensure compatibility
                        plot(axP, [x_sep, x_plat_end], [Cp(sep_idx), Cp(sep_idx)], ...
                            'Color', col, 'LineWidth', 4, 'HandleVisibility', 'off');
                        text(axP, x_sep, Cp(sep_idx), ' Plateau', 'Color', col, 'FontSize', 8, 'VerticalAlignment', 'bottom');
                    end
                    
                    % 3. LSB Patch on H Plot
                    patch(axH, [x_sep, x_rea, x_rea, x_sep], [min(ylim(axH)), min(ylim(axH)), max(ylim(axH)), max(ylim(axH))], ...
                        col, 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'HandleVisibility', 'off');
                    
                    fprintf('%s (N=%s) %s: Bubble from %.3f to %.3f\n', current_airfoil, ncrit, surfLab, x_sep, x_rea);
                end
            end
        end
    end
    legend(axH, 'Location', 'northeastoutside', 'FontSize', 7);
    linkaxes([axH, axP], 'x'); % Keep x-axis zoom synchronized
    xlim(axH, [0 1.2]);
end