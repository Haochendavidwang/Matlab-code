% flag_map_modes.m
% Professional K–L map with properties:
%  x = bending stiffness (K), y = flag length (L)
%  stability = 0 (unstable), 1 (stable)
%  mode = 0 (unstable flapping); 1 or 2 (stable flapping mode)
%
% Styling:
%  - Mode mapped to marker shape & edge color (Okabe–Ito palette)
%  - Stability mapped to fill: filled = stable, unfilled = unstable
%  - Auto axis limits with padding and degenerate-span handling
%  - Publication-friendly figure defaults

%% 1) Put your data here (each row: [K,  L,  stability,  mode])
%                    K     L    stab  mode
data = [
    0.002  2.00   0     0
    0.004  2.00   0     0
    0.006  2.00   0     0
    0.008  2.00   1     1
    0.01   2.00   1     1
    0.012  2.00   0     0
    0.014  2.00   0     0
    0.016  2.00   0     0
    0.018  2.00   0     0
    0.02   2.00   0     0

    0.002  3.00   1     1
    0.004  3.00   1     1
    0.006  3.00   1     1
    0.008  3.00   1     1
    0.01   3.00   1     1
    0.015  3.00   1     1
    0.02   3.00   1     1
    0.03   3.00   1     1
    0.04   3.00   1     1
    0.05   3.00   1     1

    0.01   4.00   0     0
    0.015  4.00   1     1
    0.02   4.00   1     1
    0.025  4.00   1     1
    0.03   4.00   1     1
    0.035  4.00   1     1
    0.04   4.00   1     1
    0.045  4.00   1     1
    0.05   4.00   1     1
    0.055  4.00   1     1
    0.06   4.00   1     1
    0.065  4.00   1     1

    0.01   5.00   0     0
    0.02   5.00   0     0
    0.03   5.00   1     2
    0.04   5.00   1     1
    0.05   5.00   1     1
    0.06   5.00   1     1
    0.08   5.00   1     1
    0.1    5.00   1     1
    0.12   5.00   1     1
    0.14   5.00   1     1
    0.16   5.00   1     1
    0.18   5.00   1     1
    0.2    5.00   1     1

    0.01   6.00   0     0
    0.01   6.00   0     0
    0.02   6.00   0     0
    0.03   6.00   1     3
    0.04   6.00   1     3
    0.05   6.00   1     3
    0.06   6.00   1     3
    0.07   6.00   1     3
    0.08   6.00   1     3
    0.09   6.00   1     3
    0.1    6.00   1     3
    0.11   6.00   1     3
    0.12   6.00   1     3
    0.13   6.00   1     3
    0.14   6.00   1     3
    0.15   6.00   1     2
    0.2    6.00   1     1
    0.25   6.00   1     1
    0.3    6.00   1     1
    0.35   6.00   1     1
    0.4    6.00   1     1
    0.45   6.00   1     1
    0.5    6.00   1     1

    0.01   7.00   0     0
    0.02   7.00   0     0
    0.03   7.00   0     0
    0.04   7.00   0     0
    0.06   7.00   0     0
    0.08   7.00   1     3
    0.1    7.00   1     3
    0.12   7.00   1     3
    0.14   7.00   1     3
    0.16   7.00   1     3
    0.16   7.00   1     3
    0.20   7.00   1     2
    0.25   7.00   1     2
    0.3    7.00   1     1
    0.35   7.00   1     1
    0.4    7.00   1     1
    0.45   7.00   1     1
    0.5    7.00   1     1
    0.55   7.00   1     1
    0.6    7.00   1     1
    0.65   7.00   1     1
    0.7    7.00   1     1
    0.75   7.00   1     1
    0.8    7.00   1     1
    0.85   7.00   1     1
    0.9    7.00   1     1
    0.95   7.00   1     1

];

% Optional: labels per row (leave {} to skip)
labels = {}; % e.g., {'A','B','C','D','E','F'}

%% 2) Appearance & layout
ms        = 70;      % marker size for scatter (points^2)
lw        = 1.2;     % marker edge width
padfrac   = 0.08;    % axis padding fraction
minspan   = 1e-6;    % degenerate-span threshold
fallback  = 1;       % span if all x or all y identical
fs_ax     = 12;      % axis font size
fs_lab    = 13;      % axis label font size
fs_title  = 14;      % title font size
fs_leg    = 11;      % legend font size

% Colorblind-safe palette (Okabe–Ito):
cBlue   = [0    114  178]/255;  % Mode 1
cOrange = [213   94    0]/255;  % Mode 2
cRed    = [213    0   50]/255;

palette = containers.Map('KeyType','double','ValueType','any');
palette(1) = cBlue;
palette(2) = cOrange;
palette(3) = cRed;

% Marker shapes per mode (extend if you later add more modes)
shapeMap = containers.Map('KeyType','double','ValueType','char');
shapeMap(1) = 'o';  % circle
shapeMap(2) = 's';  % square

%% 3) Prepare axes with automatic limits
x = data(:,1) ; y = data(:,2);
xmin = min(x); xmax = max(x); xspan = xmax - xmin;
ymin = min(y); ymax = max(y); yspan = ymax - ymin;

if xspan < minspan
    cx = (xmin + xmax)/2; xlim_lo = cx - fallback/2; xlim_hi = cx + fallback/2;
else
    xpad = padfrac * xspan; xlim_lo = xmin - xpad; xlim_hi = xmax + xpad;
end
if yspan < minspan
    cy = (ymin + ymax)/2; ylim_lo = cy - fallback/2; ylim_hi = cy + fallback/2;
else
    ypad = padfrac * yspan; ylim_lo = ymin - ypad; ylim_hi = ymax + ypad;
end

%% 4) Figure & axes
fig = figure('Color','w'); clf;
fig.Position(3:4) = [720 520]; % width x height

axes1 = axes('Parent', fig); hold(axes1, 'on'); box on; grid on; grid minor;
set(axes1, 'FontSize', fs_ax, 'LineWidth', 1.0);
xlim([xlim_lo xlim_hi]); ylim([ylim_lo ylim_hi]);
xlabel('Bending stiffness, K', 'FontSize', fs_lab);
ylabel('Flag length, L',      'FontSize', fs_lab);
title('K–L Map with Stability and Mode', 'FontSize', fs_title);

%% 5) Plot by (mode, stability)
modes = unique(data(:,4))';
stabs = [0 1]; % enforce order: 0=unstable, 1=stable

legHandles = gobjects(0);
legLabels  = strings(0);

for m = modes
    % Fallbacks if an unseen mode appears
    color = cBlue; markerShape = 'o';
    if isKey(palette, m),   color = palette(m);    end
    if isKey(shapeMap, m),  markerShape = shapeMap(m); end

    for s = stabs
        idx = data(:,4) == m & data(:,3) == s;
        if ~any(idx), continue; end

        xi = data(idx,1);
        yi = data(idx,2);

        if s == 1
            % stable: filled marker
            h = scatter(xi, yi, ms, 'Marker', markerShape, ...
                'MarkerFaceColor', color, 'MarkerEdgeColor', color, ...
                'LineWidth', lw);
            labelText = sprintf('Mode %d — Stable', m);
        else
            % unstable: hollow marker
            h = scatter(xi, yi, ms, 'Marker', markerShape, ...
                'MarkerFaceColor', 'none', 'MarkerEdgeColor', color, ...
                'LineWidth', lw);
            labelText = sprintf('Mode %d — Unstable', m);
        end

        legHandles(end+1) = h; %#ok<SAGROW>
        legLabels(end+1)  = labelText; %#ok<SAGROW>
    end
end

% Build legend (only for categories present)
if ~isempty(legHandles)
    [~, ia] = unique(legLabels, 'stable'); % avoid duplicate labels
    legend(legHandles(ia), legLabels(ia), 'Location', 'bestoutside', 'FontSize', fs_leg);
end

%% 6) Optional: per-point labels
if ~isempty(labels) && numel(labels) == size(data,1)
    dx = 0.012 * max(xlim_hi - xlim_lo, eps);
    dy = 0.012 * max(ylim_hi - ylim_lo, eps);
    for i = 1:size(data,1)
        text(x(i)+dx, y(i)+dy, labels{i}, 'FontSize', fs_ax, 'Color', [0 0 0]);
    end
end

%% 7) Optional: export (uncomment as needed)
% exportgraphics(fig, 'KL_map_modes.pdf', 'ContentType', 'vector'); % publication-quality
% exportgraphics(fig, 'KL_map_modes.png', 'Resolution', 300);

%% 8) --- NEW FIGURE: (K/L^2) vs L Map ---
% This section generates a second figure with K/L^2 on the x-axis.
% It re-uses all settings (colors, sizes) from the first plot.

% 8.1) Calculate new X-data and its limits
x_new = data(:,1) ./ (data(:,2).^2); % K/L^2
y_new = data(:,2);                   % L (same as y)

xmin_new = min(x_new); xmax_new = max(x_new); xspan_new = xmax_new - xmin_new;
% y-limits are identical to the first plot (ylim_lo, ylim_hi)

if xspan_new < minspan
    cx_new = (xmin_new + xmax_new)/2; 
    xlim_lo_new = cx_new - fallback/2; 
    xlim_hi_new = cx_new + fallback/2;
else
    xpad_new = padfrac * xspan_new; 
    xlim_lo_new = xmin_new - xpad_new; 
    xlim_hi_new = xmax_new + xpad_new;
end

% 8.2) Create new Figure & axes
fig2 = figure('Color','w'); clf;
fig2.Position(3:4) = [720 520]; % width x height

axes2 = axes('Parent', fig2); 
hold(axes2, 'on'); box on; grid on; grid minor;
set(axes2, 'FontSize', fs_ax, 'LineWidth', 1.0);
xlim(axes2, [xlim_lo_new xlim_hi_new]); 
ylim(axes2, [ylim_lo ylim_hi]); % Re-use original y-limits

% Set new labels using TeX interpreter for K/L^2
xlabel(axes2, 'Scaled Stiffness, K/L^2', 'FontSize', fs_lab, 'Interpreter', 'tex');
ylabel(axes2, 'Flag length, L',          'FontSize', fs_lab);
title(axes2,  '(K/L^2)–L Map with Stability and Mode', 'FontSize', fs_title, 'Interpreter', 'tex');

% 8.3) Plot by (mode, stability) - same loop, new x-data
% We must re-initialize legend handles for the new figure
legHandles2 = gobjects(0);
legLabels2  = strings(0);

for m = modes % 'modes' is from the first plot's setup
    % Fallbacks if an unseen mode appears
    color = cBlue; markerShape = 'o';
    if isKey(palette, m),   color = palette(m);    end
    if isKey(shapeMap, m),  markerShape = shapeMap(m); end

    for s = stabs % 'stabs' is from the first plot's setup
        idx = data(:,4) == m & data(:,3) == s;
        if ~any(idx), continue; end

        % --- Use the NEW x-data ---
        xi = x_new(idx); 
        % --- Use the original y-data ---
        yi = data(idx,2);

        if s == 1
            % stable: filled marker
            % NOTE: We explicitly plot to axes2
            h = scatter(axes2, xi, yi, ms, 'Marker', markerShape, ...
                'MarkerFaceColor', color, 'MarkerEdgeColor', color, ...
                'LineWidth', lw);
            labelText = sprintf('Mode %d — Stable', m);
        else
            % unstable: hollow marker
            % NOTE: We explicitly plot to axes2
            h = scatter(axes2, xi, yi, ms, 'Marker', markerShape, ...
                'MarkerFaceColor', 'none', 'MarkerEdgeColor', color, ...
                'LineWidth', lw);
            labelText = sprintf('Mode %d — Unstable', m);
        end

        legHandles2(end+1) = h; %#ok<SAGROW>
        legLabels2(end+1)  = labelText; %#ok<SAGROW>
    end
end

% 8.4) Build legend for the new figure
if ~isempty(legHandles2)
    [~, ia] = unique(legLabels2, 'stable'); % avoid duplicate labels
    % NOTE: We explicitly add legend to axes2
    legend(axes2, legHandles2(ia), legLabels2(ia), 'Location', 'bestoutside', 'FontSize', fs_leg);
end

% 8.5) Optional: per-point labels for new figure
if ~isempty(labels) && numel(labels) == size(data,1)
    % Calculate new dx offset based on new x-axis span
    dx_new = 0.012 * max(xlim_hi_new - xlim_lo_new, eps);
    % dy is the same as before
    dy = 0.012 * max(ylim_hi - ylim_lo, eps); 
    
    for i = 1:size(data,1)
        % NOTE: We plot text to axes2 using x_new
        text(axes2, x_new(i)+dx_new, y_new(i)+dy, labels{i}, 'FontSize', fs_ax, 'Color', [0 0 0]);
    end
end

%% 9) Optional: export NEW figure (uncomment as needed)
% exportgraphics(fig2, 'KL2_L_map_modes.pdf', 'ContentType', 'vector');
% exportgraphics(fig2, 'KL2_L_map_modes.png', 'Resolution', 300);
