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

colors = lines(size(cases,1)); 
figure('Color', 'w', 'Position', [100, 100, 1100, 700]);
hold on; grid on; box on;

% Dummy plots for legend
plot(NaN,NaN,'ko','MarkerFaceColor','k','DisplayName','Transition');
plot(NaN,NaN,'rs','MarkerFaceColor','r','DisplayName','Separation (Cf <= 0)');

fprintf('\n--- Boundary Layer Diagnostics ---\n')

%% ----------------------------------------------------
% LOOP THROUGH CASES
%% ----------------------------------------------------
for i = 1:size(cases,1)
    filename = cases{i,1};
    airfoil  = cases{i,2};
    ncrit    = cases{i,3};
    
    if ~exist(filename, 'file'), continue; end
    data = readmatrix(filename, 'FileType', 'text');
    
    x_raw  = data(:,2);
    y_raw  = data(:,3);
    Cf_raw = data(:,7);
    H_raw  = data(:,8);

    surfaces = {'upper', 'lower'};
    
    for s = 1:2
        if strcmp(surfaces{s}, 'upper')
            idx = y_raw >= 0;
            lStyle = '-'; 
            surfLab = 'Up';
        else
            idx = y_raw < 0;
            lStyle = '--'; 
            surfLab = 'Lo';
        end
        
        x  = x_raw(idx);
        Cf = Cf_raw(idx);
        H  = H_raw(idx);
        
        % Sort to prevent "zig-zag" lines
        [x, sortIdx] = sort(x);
        Cf = Cf(sortIdx);
        H = H(sortIdx);
        
        % Plot the Shape Factor H
        fullLabel = sprintf('%s (N=%s) %s', airfoil, ncrit, surfLab);
        plot(x, H, lStyle, 'LineWidth', 1.8, 'Color', colors(i,:), 'DisplayName', fullLabel);

        % ---- MARK TRANSITION ----
        tr_idx = find(H(1:end-1) > 2.2 & H(2:end) < 1.8, 1);
        if ~isempty(tr_idx)
            % Plot the marker ON the H line
            mTr = plot(x(tr_idx), H(tr_idx), 'ko', 'MarkerFaceColor', 'y', 'MarkerSize', 8, 'LineWidth', 1.2, 'HandleVisibility', 'off');
            uistack(mTr, 'top');
            fprintf('%s (N=%s) %s: Transition at x/c = %.3f\n', airfoil, ncrit, surfLab, x(tr_idx))
        else
            fprintf('%s (N=%s) %s: No transition detected\n', airfoil, ncrit, surfLab)
        end

        % ---- MARK SEPARATION ----
        % Logic: Cf drops to or below 0
        sep_idx = find(Cf <= 1e-6, 1); 
        if ~isempty(sep_idx)
            % CRITICAL FIX: Plot the marker ON the H line at the separation x-coordinate
            plot(x(sep_idx), H(sep_idx), 'rs', 'MarkerFaceColor', 'r', 'MarkerSize', 7, 'HandleVisibility', 'off');
            fprintf('%s (N=%s) %s: Separation at x/c = %.3f\n', airfoil, ncrit, surfLab, x(sep_idx))
        else
            fprintf('%s (N=%s) %s: No separation detected\n', airfoil, ncrit, surfLab)
        end

        % ---- DETECT KEY PHYSICS ----
        sep_idx = find(Cf <= 1e-6, 1); 
        tr_idx = find(H(1:end-1) > 2.2 & H(2:end) < 1.8, 1);

        % ---- HIGHLIGHT LAMINAR SEPARATION BUBBLE (LSB) ----
        if ~isempty(sep_idx)
            % Find reattachment: The first point AFTER separation where H drops 
            % to turbulent levels (H < 1.8)
            reattach_idx = find(H(sep_idx:end) < 1.8, 1) + sep_idx - 1;
            
            if ~isempty(reattach_idx)
                x_start = x(sep_idx);
                x_end   = x(reattach_idx);
                
                % Define the pale patch area
                x_patch = [x_start, x_end, x_end, x_start];
                y_patch = [1, 1, 5, 5]; % Covers the H-axis range
                
                % Use a more visible patch style
                patch(x_patch, y_patch, colors(i,:), ...
                    'FaceAlpha', 0.2, ...           % Slightly darker for visibility
                    'EdgeColor', colors(i,:), ...   % Solid edge to define boundaries
                    'LineWidth', 1, ...
                    'LineStyle', '--', ...
                    'HandleVisibility', 'off');
                
                % Mark transition at the reattachment point for clarity
                plot(x(reattach_idx), H(reattach_idx), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'HandleVisibility', 'off');
                
                % Add a text label for length
                len = (x_end - x_start) * 100;
                text(x_start, 4.8, sprintf('%.1f%%', len), 'Color', colors(i,:), ...
                    'FontSize', 8, 'FontWeight', 'bold');
            end
        end
    end
end

%% ----------------------------------------------------
% FINAL PLOT SETTINGS
%% ----------------------------------------------------
xlabel('Chordwise Position (x/c)', 'FontSize', 12, 'FontWeight', 'bold'); 
ylabel('Shape Factor (H)', 'FontSize', 12, 'FontWeight', 'bold');
title('Shape Factor Comparison with Transition and Separation Markers');
legend('Location', 'northeastoutside', 'FontSize',8);
xlim([0 1.5]); % Focus on airfoil; ignore wake
ylim([1 5]);    % Standard range for H (1.4=turbulent, 2.6=laminar, >3.5=separated)

