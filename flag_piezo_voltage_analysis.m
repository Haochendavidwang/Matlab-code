%% MATLAB代码：计算多个旗帜的压电电压（曲率积分）
% 仅考虑pinned支撑方式
% 计算每个旗帜的曲率积分（代表压电输出电压）并进行比较

clear all; close all; clc;

%% 1. 配置
% =========================================================================
% 定义要分析的文件配置
% - 'filename':    数据文件的完整路径
% - 'L':           旗帜长度参数
% - 'DisplayName': 用于图例的名称
% - 't_min':       开始时间
% - 't_max':       结束时间
% =========================================================================
configs = [
    struct('filename', 'F:\code fortran\cylinderstring\D1L2\D1L2_0.01\SXYM.dat', ...
           'L', 2, 'DisplayName', 'L=2, mode=0.01', 't_min', 250, 't_max', 300);
           
    struct('filename', 'F:\code fortran\cylinderstring\D1L3\D1L3_0.02\SXYM.dat', ...
           'L', 3, 'DisplayName', 'L=3, mode=0.02', 't_min', 250, 't_max', 300);
           
    struct('filename', 'F:\code fortran\cylinderstring\D1L4\D1L4_0.04\SXYM.dat', ...
           'L', 4, 'DisplayName', 'L=4, mode=0.04', 't_min', 250, 't_max', 300); 
          
    struct('filename', 'F:\code fortran\cylinderstring\D1L5\D1L5_0.03\SXYM.dat', ...
           'L', 5, 'DisplayName', 'L=5, mode2 (0.03)', 't_min', 250, 't_max', 300);

    struct('filename', 'F:\code fortran\cylinderstring\D1L5\D1L5_0.12\SXYM.dat', ...
           'L', 5, 'DisplayName', 'L=5, mode1 (0.12)', 't_min', 150, 't_max', 200);
           
    struct('filename', 'F:\code fortran\cylinderstring\D1L6\D1L6_0.05\SXYM.dat', ...
           'L', 6, 'DisplayName', 'L=6, mode2 (0.05)', 't_min', 150, 't_max', 200);

    struct('filename', 'F:\code fortran\cylinderstring\D1L6\D1L6_0.4\SXYM.dat', ...
           'L', 6, 'DisplayName', 'L=6, mode1 (0.4)', 't_min', 150, 't_max', 200);

    struct('filename', 'F:\code fortran\cylinderstring\D1L7\D1L7_0.08\SXYM.dat', ...
           'L', 7, 'DisplayName', 'L=7, mode2 (0.08)', 't_min', 150, 't_max', 200);
];

% 分析参数
track_point_index = -1;  % -1表示末尾点，用于轨迹追踪

%% 2. 数据存储初始化
N_files = length(configs);
all_results = cell(N_files, 1);

%% 3. 主处理循环
fprintf('开始处理 %d 个文件...\n', N_files);
fprintf('========================================\n');

for i = 1:N_files
    % 获取当前文件配置
    full_path = configs(i).filename;
    L = configs(i).L;
    t_min_local = configs(i).t_min;
    t_max_local = configs(i).t_max;
    display_name = configs(i).DisplayName;
    
    % 计算期望的点数（pinned情况）
    expected_points = L * 32 + 1;
    
    fprintf('正在处理文件 %d/%d: %s\n', i, N_files, display_name);
    fprintf('  路径: %s\n', full_path);
    fprintf('  L=%d, 期望点数=%d, 时间窗口=[%.2f, %.2f]\n', ...
            L, expected_points, t_min_local, t_max_local);
    
    % 读取旗帜数据（第二个ZONE）
    [time_points, X_tracked, Y_tracked, allX, allY] = ...
        readFlagData(full_path, expected_points, track_point_index);
    
    if isempty(time_points)
        warning('文件 %s 未能读取到有效数据', full_path);
        continue;
    end
    
    % 应用时间窗口
    mask = (time_points >= t_min_local) & (time_points <= t_max_local);
    
    if sum(mask) == 0
        warning('时间窗口 [%.2f, %.2f] 内没有数据', t_min_local, t_max_local);
        continue;
    end
    
    % 计算曲率积分
    fprintf('  计算曲率积分...\n');
    [int_curv_sq, int_curv, local_curv] = ...
        calculateCurvatureIntegral(allX, allY, time_points, mask);
    
    % 存储结果
    results = struct();
    results.config = configs(i);
    results.time = time_points(mask);
    results.X_end = X_tracked(mask);
    results.Y_end = Y_tracked(mask);
    results.allX = allX(mask);
    results.allY = allY(mask);
    results.int_curv_sq = int_curv_sq;  % 弯曲能量
    results.int_curv = int_curv;        % 压电电压（考虑正负抵消）
    results.local_curv = local_curv;
    
    all_results{i} = results;
    
    % 计算统计信息
    max_voltage = max(abs(int_curv));
    rms_voltage = sqrt(mean(int_curv.^2));
    
    fprintf('  完成！有效时间点: %d\n', sum(mask));
    fprintf('  压电电压统计:\n');
    fprintf('    最大值: %.6f\n', max_voltage);
    fprintf('    RMS值: %.6f\n', rms_voltage);
    fprintf('----------------------------------------\n');
end

fprintf('所有文件处理完毕！\n\n');

%% 4. 绘制压电电压时程对比图
figure('Name', '旗帜压电电压时程对比', 'Position', [100, 100, 1400, 800]);

% 子图1：压电电压（∫κds）
subplot(2, 1, 1);
hold on;
colors = lines(N_files);
legend_entries = {};

for i = 1:N_files
    if ~isempty(all_results{i})
        results = all_results{i};
        plot(results.time, results.int_curv, '-', ...
             'LineWidth', 1.5, 'Color', colors(i,:), ...
             'DisplayName', results.config.DisplayName);
        legend_entries{end+1} = results.config.DisplayName;
    end
end

xlabel('时间 t', 'FontSize', 12);
ylabel('压电输出电压 ∫κds', 'FontSize', 12);
title('旗帜压电输出电压对比', 'FontSize', 14, 'FontWeight', 'bold');
legend(legend_entries, 'Location', 'bestoutside', 'FontSize', 10);
grid on;
box on;

%% 5. 统计分析和性能对比
fprintf('\n========== 压电电压性能统计分析 ==========\n');
fprintf('%-20s | %-10s | %-10s | %-10s | %-10s\n', ...
        '配置', '平均|电压|', '最大|电压|', 'RMS电压');
fprintf('---------------------------------------------------------------------\n');

% 收集统计数据用于排序
stats_data = [];
valid_configs = {};

for i = 1:N_files
    if ~isempty(all_results{i})
        results = all_results{i};
     
        max_abs_voltage = max(abs(results.int_curv));
        rms_voltage = sqrt(mean(results.int_curv.^2));

        
        fprintf('%-20s | %10.6f | %10.6f \n', ...
                results.config.DisplayName, ...
                max_abs_voltage, rms_voltage);
        
        stats_data(end+1, :) = [max_abs_voltage, rms_voltage];
        valid_configs{end+1} = results.config.DisplayName;
    end
end

%% 6. 创建性能对比柱状图
if ~isempty(stats_data)
    figure('Name', '压电电压性能对比', 'Position', [150, 150, 1200, 600]);
    
    % 创建柱状图
    subplot(1, 2, 1);
    bar_data = [stats_data(:,1), stats_data(:,2)];
    b = bar(bar_data);
    
    % 设置颜色

    b(1).FaceColor = [0.8, 0.2, 0.2];  % 最大值 - 红色
    b(2).FaceColor = [0.2, 0.7, 0.2];  % RMS值 - 绿色
    
    set(gca, 'XTickLabel', valid_configs);
    xtickangle(45);
    ylabel('电压幅度', 'FontSize', 12);
    title('压电电压性能指标对比', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'最大|电压|', 'RMS电压'}, 'Location', 'best');
    grid on;
    
    % 创建雷达图（如果配置数少于等于8个）
    if length(valid_configs) <= 8
        subplot(1, 2, 2);
        
        % 归一化数据（0-1范围）
        norm_data = stats_data;
        for j = 1:size(norm_data, 2)
            max_val = max(norm_data(:, j));
            if max_val > 0
                norm_data(:, j) = norm_data(:, j) / max_val;
            end
        end
        
        % 绘制雷达图
        theta = linspace(0, 2*pi, size(norm_data, 1) + 1);
        
        hold on;
        metrics = {'平均电压', '最大电压', 'RMS电压', '振荡幅度'};
        colors_radar = [0.2, 0.4, 0.8; 0.8, 0.2, 0.2; 0.2, 0.7, 0.2; 0.7, 0.5, 0.1];
        
        for j = 1:4
            r = [norm_data(:, j); norm_data(1, j)]';
            polarplot(theta, r, 'o-', 'LineWidth', 2, ...
                     'Color', colors_radar(j,:), ...
                     'MarkerFaceColor', colors_radar(j,:));
        end
        
        % 设置极坐标轴
        ax = gca;
        ax.ThetaTick = (0:length(valid_configs)-1) * 360 / length(valid_configs);
        ax.ThetaTickLabel = valid_configs;
        ax.RLim = [0, 1];
        title('归一化性能雷达图', 'FontSize', 14, 'FontWeight', 'bold');
        legend(metrics, 'Location', 'bestoutside');
    end
end


%% 辅助函数

% 读取旗帜数据（第二个ZONE）
function [time_points, X_tracked, Y_tracked, allX, allY] = readFlagData(filename, expected_points, track_index)
    % 读取文件
    fid = fopen(filename, 'r');
    if fid == -1
        warning('无法打开文件: %s', filename);
        time_points = [];
        X_tracked = [];
        Y_tracked = [];
        allX = [];
        allY = [];
        return;
    end
    
    try
        file_content = textscan(fid, '%s', 'Delimiter', '\n');
        lines = file_content{1};
        fclose(fid);
    catch ME
        fclose(fid);
        warning('读取文件时发生错误: %s', ME.message);
        time_points = [];
        X_tracked = [];
        Y_tracked = [];
        allX = [];
        allY = [];
        return;
    end
    
    % 查找ZONE
    zone_starts = find(contains(lines, 'ZONE T='));
    n_zones = length(zone_starts);
    
    % 每两个ZONE为一组
    n_time_steps = floor(n_zones / 2);
    
    % 预分配
    time_points = zeros(n_time_steps, 1);
    X_tracked = zeros(n_time_steps, 1);
    Y_tracked = zeros(n_time_steps, 1);
    allX = cell(n_time_steps, 1);
    allY = cell(n_time_steps, 1);
    
    valid_steps = 0;
    
    % 处理每个时间步（读取第二个ZONE - 旗帜数据）
    for i = 1:n_time_steps
        flag_zone_idx = 2 * i;  % 第二个ZONE
        
        if flag_zone_idx > n_zones
            break;
        end
        
        % 提取时间
        zone_line = lines{zone_starts(flag_zone_idx)};
        time_match = regexp(zone_line, 'T=".*_\s*([\d\.\+\-E]+)"', 'tokens');
        
        if isempty(time_match)
            continue;
        end
        
        current_time = str2double(time_match{1}{1});
        
        % 找到数据开始位置
        data_start = zone_starts(flag_zone_idx);
        while data_start <= length(lines) && ~contains(lines{data_start}, 'I=')
            data_start = data_start + 1;
        end
        data_start = data_start + 1;
        
        % 确定数据结束位置
        if flag_zone_idx < n_zones
            data_end = zone_starts(flag_zone_idx + 1) - 1;
        else
            data_end = length(lines);
        end
        
        % 读取数据
        [X_data, Y_data] = readZoneData(lines, data_start, data_end, expected_points);
        
        % 验证数据
        np = min(length(X_data), length(Y_data));
        if np < expected_points
            continue;
        end
        
        % 存储数据
        valid_steps = valid_steps + 1;
        time_points(valid_steps) = current_time;
        
        allX{valid_steps} = X_data(1:expected_points);
        allY{valid_steps} = Y_data(1:expected_points);
        
        % 追踪特定点
        if track_index == -1
            track_idx = expected_points;
        elseif track_index > 0 && track_index <= expected_points
            track_idx = track_index;
        else
            track_idx = expected_points;
        end
        
        X_tracked(valid_steps) = X_data(track_idx);
        Y_tracked(valid_steps) = Y_data(track_idx);
    end
    
    % 截取有效数据
    time_points = time_points(1:valid_steps);
    X_tracked = X_tracked(1:valid_steps);
    Y_tracked = Y_tracked(1:valid_steps);
    allX = allX(1:valid_steps);
    allY = allY(1:valid_steps);
end

% 读取ZONE数据
function [X_data, Y_data] = readZoneData(lines, data_start, data_end, expected_points)
    X_data = [];
    Y_data = [];
    
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
    
    % 分割数据
    if length(all_values) >= 2*expected_points
        X_data = all_values(1:expected_points);
        Y_data = all_values(expected_points+1:2*expected_points);
    elseif length(all_values) >= expected_points
        mid_point = floor(length(all_values)/2);
        X_data = all_values(1:mid_point);
        Y_data = all_values(mid_point+1:min(2*mid_point, end));
    end
    
    X_data = X_data(:);
    Y_data = Y_data(:);
end

% 计算曲率积分
function [int_curv_sq, int_curv, local_curv] = calculateCurvatureIntegral(allX, allY, time, mask)
    indices = find(mask);
    n_times = length(indices);
    int_curv_sq = zeros(n_times, 1);
    int_curv = zeros(n_times, 1);
    local_curv = cell(n_times, 1);
    
    for i = 1:n_times
        idx = indices(i);
        x = allX{idx};
        y = allY{idx};
        
        n_pts = length(x);
        if n_pts >= 3
            % 计算弧长参数
            ds = sqrt(diff(x).^2 + diff(y).^2);
            s = [0; cumsum(ds)];
            
            % 初始化曲率数组
            kappa = zeros(n_pts, 1);
            
            % 内部点使用中心差分
            for j = 2:n_pts-1
                x_pts = x(j-1:j+1);
                y_pts = y(j-1:j+1);
                s_pts = s(j-1:j+1);
                
                if length(unique(s_pts)) == 3
                    % 二阶多项式拟合
                    px = polyfit(s_pts, x_pts, 2);
                    py = polyfit(s_pts, y_pts, 2);
                    
                    % 计算导数
                    dx_ds = 2*px(1)*s(j) + px(2);
                    dy_ds = 2*py(1)*s(j) + py(2);
                    d2x_ds2 = 2*px(1);
                    d2y_ds2 = 2*py(1);
                    
                    % 计算曲率
                    numerator = dx_ds * d2y_ds2 - dy_ds * d2x_ds2;
                    denominator = (dx_ds^2 + dy_ds^2)^(3/2);
                    
                    if denominator > 1e-10
                        kappa(j) = numerator / denominator;
                    end
                end
            end
            
            % 边界点处理
            % 第一个点
            if n_pts >= 3
                x_pts = x(1:3);
                y_pts = y(1:3);
                s_pts = s(1:3);
                if length(unique(s_pts)) == 3
                    px = polyfit(s_pts, x_pts, 2);
                    py = polyfit(s_pts, y_pts, 2);
                    dx_ds = 2*px(1)*s(1) + px(2);
                    dy_ds = 2*py(1)*s(1) + py(2);
                    d2x_ds2 = 2*px(1);
                    d2y_ds2 = 2*py(1);
                    numerator = dx_ds * d2y_ds2 - dy_ds * d2x_ds2;
                    denominator = (dx_ds^2 + dy_ds^2)^(3/2);
                    if denominator > 1e-10
                        kappa(1) = numerator / denominator;
                    end
                end
            end
            
            % 最后一个点
            if n_pts >= 3
                x_pts = x(end-2:end);
                y_pts = y(end-2:end);
                s_pts = s(end-2:end);
                if length(unique(s_pts)) == 3
                    px = polyfit(s_pts, x_pts, 2);
                    py = polyfit(s_pts, y_pts, 2);
                    dx_ds = 2*px(1)*s(end) + px(2);
                    dy_ds = 2*py(1)*s(end) + py(2);
                    d2x_ds2 = 2*px(1);
                    d2y_ds2 = 2*py(1);
                    numerator = dx_ds * d2y_ds2 - dy_ds * d2x_ds2;
                    denominator = (dx_ds^2 + dy_ds^2)^(3/2);
                    if denominator > 1e-10
                        kappa(end) = numerator / denominator;
                    end
                end
            end
            
            % 计算积分
            kappa_sq = kappa.^2;
            int_curv_sq(i) = trapz(s, kappa_sq);  % 弯曲能量
            int_curv(i) = trapz(s, kappa);        % 压电电压（考虑正负）
            
            local_curv{i} = struct('s', s, 'kappa', kappa);
        end
    end
end

fprintf('\n分析完成！请查看生成的图表。\n');
