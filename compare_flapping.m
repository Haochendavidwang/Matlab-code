%% MATLAB代码：比较多个DAT文件中的旗帜轮廓 (V3)
%
% 关键改动：
% 1. 'configs' 结构现在包含 't_min' 和 't_max'，允许为每个文件
%    指定 *独立的时间窗口*。
% 2. 移除全局时间窗口设置。
% 3. 数据处理循环现在存储每个时间步的 *完整X, Y轮廓*。
% 4. 绘图部分被修改为在每个子图中绘制轮廓演化 (X-Y 叠加图)。
% 5. 使用 'axis equal' 和 'linkaxes(..., 'xy')' 来统一所有子图的
%    空间坐标系，以便于比较。

% 清理工作空间
clear all; close all; clc;

%% 1. 配置
% =========================================================================
% 在这里定义你要比较的文件
%
% - 'filename':    数据文件的 *完整路径*
% - 'L':           旗帜长度参数
% - 'DisplayName': 用于图例的名称
% - 't_min':       此文件的开始时间
% - 't_max':       此文件的结束时间
% =========================================================================
configs = [
    struct('filename', 'F:\code fortran\cylinderstring\D1L2\D1L2_0.01\SXYM.dat', ...
           'L', 2, 'DisplayName', 'L=2,0.01', 't_min', 270, 't_max', 300, 'mode', 1);
           
    struct('filename', 'F:\code fortran\cylinderstring\D1L3\D1L3_0.02\SXYM.dat', ...
           'L', 3, 'DisplayName', 'L=3,0.02', 't_min', 270, 't_max', 300, 'mode', 1);
           
    struct('filename', 'F:\code fortran\cylinderstring\D1L4\D1L4_0.04\SXYM.dat', ...
           'L', 4, 'DisplayName', 'L=4,0.04', 't_min', 270, 't_max', 300, 'mode', 1); 
          
    struct('filename', 'F:\code fortran\cylinderstring\D1L5\D1L5_0.03\SXYM.dat', ...
           'L', 5, 'DisplayName', 'L=5,0.03', 't_min', 270, 't_max', 300, 'mode', 2);

    struct('filename', 'F:\code fortran\cylinderstring\D1L5\D1L5_0.12\SXYM.dat', ...
           'L', 5, 'DisplayName', 'L=5,0.12', 't_min', 170, 't_max', 200, 'mode', 1);
           
    struct('filename', 'F:\code fortran\cylinderstring\D1L6\D1L6_0.05\SXYM.dat', ...
           'L', 6, 'DisplayName', 'L=6,0.05', 't_min', 170, 't_max', 200, 'mode', 3);

    struct('filename', 'F:\code fortran\cylinderstring\D1L6\D1L6_0.15\SXYM.dat', ...
           'L', 6, 'DisplayName', 'L=6,0.4', 't_min', 170, 't_max', 200, 'mode', 2);

    struct('filename', 'F:\code fortran\cylinderstring\D1L6\D1L6_0.4\SXYM.dat', ...
           'L', 6, 'DisplayName', 'L=6,0.15', 't_min', 270, 't_max', 300, 'mode', 1);

    struct('filename', 'F:\code fortran\cylinderstring\D1L7\D1L7_0.08\SXYM.dat', ...
           'L', 7, 'DisplayName', 'L=7,0.08', 't_min', 170, 't_max', 200, 'mode', 3);

    struct('filename', 'F:\code fortran\cylinderstring\D1L7\D1L7_0.2\SXYM.dat', ...
           'L', 7, 'DisplayName', 'L=7,0.2', 't_min', 270, 't_max', 300, 'mode', 2);

    struct('filename', 'F:\code fortran\cylinderstring\D1L7\D1L7_0.7\SXYM.dat', ...
           'L', 7, 'DisplayName', 'L=7,0.7', 't_min', 270, 't_max', 300, 'mode', 1);
];

% 绘图参数
max_display_lines = 100;    % 每个子图最大显示轮廓线数量
profile_line_width = 0.5;   % 轮廓线宽 (设小一点以便看清)


%% 2. 数据存储
N_files = length(configs);
all_results = cell(N_files, 1);  % 使用元胞数组存储每个文件的结果

%% 3. 主处理循环 (循环遍历每个文件)
fprintf('开始处理 %d 个文件...\n', N_files);
fprintf('========================================\n');

for i = 1:N_files
    % --- 获取当前文件的配置 ---
    full_path = configs(i).filename;
    L = configs(i).L;
    t_min_local = configs(i).t_min;
    t_max_local = configs(i).t_max;
    
    % --- 关键逻辑：根据L计算期望的点数 ---
    expected_points = L * 32 + 1;
    
    fprintf('正在处理文件: %s\n(L=%d, expected_points=%d, T=[%.2f, %.2f])\n', ...
            full_path, L, expected_points, t_min_local, t_max_local);

    % --- 1. 读取文件 ---
    fid = fopen(full_path, 'r');
    if fid == -1
        warning('无法打开文件 %s，跳过此文件', full_path);
        continue;
    end
    try
        file_content = textscan(fid, '%s', 'Delimiter', '\n');
        lines = file_content{1};
        fclose(fid);
    catch ME
        fclose(fid);
        warning('读取文件 %s 时发生错误: %s', full_path, ME.message);
        continue;
    end

    % --- 2. 查找ZONE并预分配 ---
    zone_starts = find(contains(lines, 'ZONE T='));
    n_zones = length(zone_starts);

    valid_time_zones = 0;
    for k = 1:2:n_zones  % 每两个ZONE为一组
        if k+1 <= n_zones
            valid_time_zones = valid_time_zones + 1;
        end
    end
    
    if valid_time_zones == 0
        warning('文件 %s 中未找到有效的ZONE对。', full_path);
        continue;
    end

    % 预分配数组 (现在存储完整的轮廓)
    time_points = zeros(valid_time_zones, 1);
    allX_local = cell(valid_time_zones, 1);
    allY_local = cell(valid_time_zones, 1);
    valid_zones = 0;

    % --- 3. 处理每组ZONE (假定总是读取第二组ZONE) ---
    zone_idx = 1;
    while zone_idx <= n_zones
        if zone_idx+1 > n_zones
            break;
        end
        target_zone_idx = zone_idx + 1;
        
        zone_line = lines{zone_starts(target_zone_idx)};
        time_match = regexp(zone_line, 'T=".*_\s*([\d\.\+\-E]+)"', 'tokens');
        
        if isempty(time_match)
            fprintf('警告 (文件 %s): 无法解析第%d个ZONE的时间\n', full_path, target_zone_idx);
            zone_idx = zone_idx + 2;
            continue;
        end
        
        current_time = str2double(time_match{1}{1});
        
        data_start = zone_starts(target_zone_idx);
        while data_start <= length(lines) && ~contains(lines{data_start}, 'I=')
            data_start = data_start + 1;
        end
        data_start = data_start + 1;
        
        if target_zone_idx < n_zones
            data_end = zone_starts(target_zone_idx+1) - 1;
        else
            data_end = length(lines);
        end
        
        [X_data, Y_data] = readZoneDataNew(lines, data_start, data_end, expected_points);
        
        np = min(length(X_data), length(Y_data));
        if np < 2
            fprintf('警告 (文件 %s): 时间 %.6f 数据点不足\n', full_path, current_time);
            zone_idx = zone_idx + 2;
            continue;
        end
        
        % 存储有效数据 (完整的轮廓)
        valid_zones = valid_zones + 1;
        time_points(valid_zones) = current_time;
        allX_local{valid_zones} = X_data(1:np);
        allY_local{valid_zones} = Y_data(1:np);
        
        zone_idx = zone_idx + 2;
    end

    % --- 4. 截取并存储结果 ---
    time_points = time_points(1:valid_zones);
    allX_local = allX_local(1:valid_zones);
    allY_local = allY_local(1:valid_zones);

    % 应用 *局部的* 时间窗口
    mask = (time_points >= t_min_local) & (time_points <= t_max_local);

    % 存储最终结果
    results.time = time_points(mask);
    results.allX = allX_local(mask);
    results.allY = allY_local(mask);
    results.config = configs(i);
    
    all_results{i} = results;
    fprintf('文件 %s 处理完成, 找到 %d 个有效时间点。\n', full_path, sum(mask));
    fprintf('----------------------------------------\n');
end

fprintf('所有文件处理完毕！\n');

%% 4. 绘图部分 (N x 1 轮廓子图)
fprintf('正在生成比较图...\n');

valid_plot_count = 0;
for i = 1:N_files
    if ~isempty(all_results{i}) && ~isempty(all_results{i}.time)
        valid_plot_count = valid_plot_count + 1;
    end
end

if valid_plot_count == 0
    error('没有成功处理任何数据，无法绘图。');
end

% Create figure with tiledlayout for no gaps
hFig = figure('Name', '旗帜轮廓演化比较', 'Position', [100, 50, 900, 150 * valid_plot_count]);

% Create tiledlayout with no padding
t = tiledlayout(valid_plot_count, 1, 'TileSpacing', 'none', 'Padding', 'compact');

ax = gobjects(valid_plot_count, 1);

plot_idx = 0;
for i = 1:N_files
    results = all_results{i};
    
    if isempty(results) || isempty(results.time)
        continue;
    end
    
    plot_idx = plot_idx + 1;

    % Create tile
    ax(plot_idx) = nexttile;
    hold(ax(plot_idx), 'on');
    
    local_allX = results.allX;
    local_allY = results.allY;
    
    % --- 选择要显示的轮廓 ---
    num_profiles = length(local_allX);
    if num_profiles == 0
        continue;
    end
    
    if num_profiles > max_display_lines
        idxs_to_plot = round(linspace(1, num_profiles, max_display_lines));
    else
        idxs_to_plot = 1:num_profiles;
    end
    
    % --- 绘制轮廓 ---
    for k = idxs_to_plot
        xk = local_allX{k};
        yk = local_allY{k};
        if length(xk) >= 2 && length(yk) >= 2
            % 全部绘制为黑色
            plot(ax(plot_idx), xk, yk, 'k-', 'LineWidth', profile_line_width);
        end
    end
    
    hold(ax(plot_idx), 'off');
    
    % --- 设置坐标轴 ---
    ylabel(['Y (', results.config.DisplayName, ')'], 'FontSize', 10, 'FontWeight', 'bold');
    grid(ax(plot_idx), 'on');
    box(ax(plot_idx), 'on');
    
    % *** 关键：使用 'axis equal' 来正确显示形状 ***
    axis(ax(plot_idx), 'equal'); 
    
    if plot_idx < valid_plot_count
        set(ax(plot_idx), 'XTickLabel', []);
    end
end

if ~isempty(ax)
    % --- 关键：链接所有X和Y轴，以便比较 ---
    linkaxes(ax, 'xy');
    
    % 为底部的图表添加X轴标签
    xlabel(ax(end), 'X 坐标', 'FontSize', 12, 'FontWeight', 'bold');
    
    % 为整个图表添加一个总标题
    sgtitle('旗帜轮廓演化比较', 'FontSize', 14, 'FontWeight', 'bold');
end

% Set the X-axis range for all linked plots
xlim(ax(1), [1, 10]);  % 示例: 设置X轴范围从 -2 到 12

fprintf('绘图完成！\n');


%% 5. 辅助函数 (与原脚本相同)
function [X_data, Y_data] = readZoneDataNew(lines, data_start, data_end, expected_points)
    % 初始化
    X_data = [];
    Y_data = [];
    
    % 读取所有数值数据
    all_values = [];
    for j = data_start:data_end
        line = strtrim(lines{j});
        
        if isempty(line) || contains(line, 'ZONE')
            break;
        end
        
        try
            values = str2double(strsplit(line));
            values = values(~isnan(values));
            if ~isempty(values)
                all_values = [all_values, values];
            end
        catch
            continue;
        end
    end
    
    % 根据期望的点数分割数据
    if length(all_values) >= 2*expected_points
        X_data = all_values(1:expected_points);
        Y_data = all_values(expected_points+1:2*expected_points);
    elseif length(all_values) >= expected_points
        mid_point = floor(length(all_values)/2);
        X_data = all_values(1:mid_point);
        Y_data = all_values(mid_point+1:min(2*mid_point, end));
    else
        fprintf('警告：数据点数量不足，期望%d对坐标，实际只有%d个值\n', ...
                expected_points, length(all_values));
    end
    
    X_data = X_data(:);
    
    % [!!] 已修正您代码中的拼写错误
    % 原始: Y_data = X_data(:); (这会使Y等于X)
    % 修正: Y_data = Y_data(:);
    Y_data = Y_data(:);
end