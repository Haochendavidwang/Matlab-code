%% MATLAB代码：计算多个旗帜的压电电压（曲率积分除以长度）
% 仅考虑pinned支撑方式
% 计算每个旗帜的曲率积分并除以旗帜长度进行归一化

clear all; close all; clc;

%% 1. 配置
configs = [
    struct('filename', 'F:\code fortran\cylinderstring\D1L2\D1L2_0.01\SXYM.dat', ...
           'L', 2, 'DisplayName', 'L=2', 't_min', 250, 't_max', 300);
           
    struct('filename', 'F:\code fortran\cylinderstring\D1L3\D1L3_0.02\SXYM.dat', ...
           'L', 3, 'DisplayName', 'L=3', 't_min', 250, 't_max', 300);
           
    struct('filename', 'F:\code fortran\cylinderstring\D1L4\D1L4_0.04\SXYM.dat', ...
           'L', 4, 'DisplayName', 'L=4', 't_min', 250, 't_max', 300); 
          
    struct('filename', 'F:\code fortran\cylinderstring\D1L5\D1L5_0.03\SXYM.dat', ...
           'L', 5, 'DisplayName', 'L=5 mode2', 't_min', 250, 't_max', 300);

    struct('filename', 'F:\code fortran\cylinderstring\D1L5\D1L5_0.12\SXYM.dat', ...
           'L', 5, 'DisplayName', 'L=5 mode1', 't_min', 150, 't_max', 200);
           
    struct('filename', 'F:\code fortran\cylinderstring\D1L6\D1L6_0.05\SXYM.dat', ...
           'L', 6, 'DisplayName', 'L=6 mode2', 't_min', 150, 't_max', 200);

    struct('filename', 'F:\code fortran\cylinderstring\D1L6\D1L6_0.4\SXYM.dat', ...
           'L', 6, 'DisplayName', 'L=6 mode1', 't_min', 150, 't_max', 200);

    struct('filename', 'F:\code fortran\cylinderstring\D1L7\D1L7_0.08\SXYM.dat', ...
           'L', 7, 'DisplayName', 'L=7 mode2', 't_min', 150, 't_max', 200);
];

%% 2. 数据处理
N_files = length(configs);
all_results = cell(N_files, 1);

fprintf('开始处理 %d 个文件...\n', N_files);
fprintf('========================================\n');

for i = 1:N_files
    % 获取配置
    full_path = configs(i).filename;
    L = configs(i).L;
    t_min_local = configs(i).t_min;
    t_max_local = configs(i).t_max;
    
    % 计算期望的点数
    expected_points = L * 32 + 1;
    
    fprintf('处理 %s (L=%d)\n', configs(i).DisplayName, L);
    
    % 读取数据
    [time_points, allX, allY] = readFlagData(full_path, expected_points);
    
    if isempty(time_points)
        warning('文件 %s 未能读取到有效数据', full_path);
        continue;
    end
    
    % 应用时间窗口
    mask = (time_points >= t_min_local) & (time_points <= t_max_local);
    
    if sum(mask) == 0
        warning('时间窗口内没有数据');
        continue;
    end
    
    % 计算曲率积分
    int_curv = calculateCurvature(allX, allY, time_points, mask);
    
    % 归一化：除以旗帜长度L
    int_curv_normalized = int_curv / L;
    
    % 存储结果
    results = struct();
    results.config = configs(i);
    results.time = time_points(mask);
    results.int_curv = int_curv;  % 原始值
    results.int_curv_normalized = int_curv_normalized;  % 归一化值
    results.L = L;
    
    all_results{i} = results;
    
    % 输出统计
    max_voltage = max(abs(int_curv_normalized));
    rms_voltage = sqrt(mean(int_curv_normalized.^2));
    
    fprintf('  最大归一化电压: %.6f\n', max_voltage);
    fprintf('  RMS归一化电压: %.6f\n', rms_voltage);
    fprintf('----------------------------------------\n');
end

fprintf('所有文件处理完毕！\n\n');

%% 3. 绘制归一化压电电压时程图
figure('Name', '归一化压电电压对比', 'Position', [100, 100, 1200, 600]);

hold on;
colors = lines(N_files);

for i = 1:N_files
    if ~isempty(all_results{i})
        results = all_results{i};
        plot(results.time, results.int_curv_normalized, '-', ...
             'LineWidth', 1.5, 'Color', colors(i,:), ...
             'DisplayName', results.config.DisplayName);
    end
end

xlabel('时间 t', 'FontSize', 12);
ylabel('归一化压电电压 (∫κds)/L', 'FontSize', 12);
title('旗帜归一化压电输出电压对比', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'bestoutside', 'FontSize', 10);
grid on;
box on;

%% 4. 统计分析表格
fprintf('\n========== 归一化压电电压统计 ==========\n');
fprintf('%-15s | %-8s | %-12s | %-12s\n', ...
        '配置', '长度L', '最大|电压|/L', 'RMS电压/L');
fprintf('--------------------------------------------------\n');

stats_data = [];
valid_configs = {};

for i = 1:N_files
    if ~isempty(all_results{i})
        results = all_results{i};
        
        max_abs_voltage = max(abs(results.int_curv_normalized));
        rms_voltage = sqrt(mean(results.int_curv_normalized.^2));
        
        fprintf('%-15s | %8d | %12.6f | %12.6f\n', ...
                results.config.DisplayName, ...
                results.L, ...
                max_abs_voltage, ...
                rms_voltage);
        
        stats_data(end+1, :) = [max_abs_voltage, rms_voltage];
        valid_configs{end+1} = results.config.DisplayName;
    end
end

%% 5. 性能对比柱状图
if ~isempty(stats_data)
    figure('Name', '归一化压电电压性能对比', 'Position', [150, 150, 800, 600]);
    
    bar_data = stats_data;
    b = bar(bar_data);
    
    % 设置颜色
    b(1).FaceColor = [0.8, 0.2, 0.2];  % 最大值 - 红色
    b(2).FaceColor = [0.2, 0.7, 0.2];  % RMS值 - 绿色
    
    set(gca, 'XTickLabel', valid_configs);
    xtickangle(45);
    ylabel('归一化电压幅度', 'FontSize', 12);
    title('归一化压电电压性能指标对比', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'最大|电压|/L', 'RMS电压/L'}, 'Location', 'best');
    grid on;
end

%% 6. 最优配置分析
fprintf('\n========== 最优配置分析 ==========\n');

if ~isempty(stats_data)
    [max_val, idx_max] = max(stats_data(:, 1));
    [rms_val, idx_rms] = max(stats_data(:, 2));
    
    fprintf('最高归一化峰值电压: %s (%.6f)\n', valid_configs{idx_max}, max_val);
    fprintf('最高归一化RMS电压: %s (%.6f)\n', valid_configs{idx_rms}, rms_val);
end

%% 辅助函数

% 简化的读取函数
function [time_points, allX, allY] = readFlagData(filename, expected_points)
    fid = fopen(filename, 'r');
    if fid == -1
        time_points = [];
        allX = [];
        allY = [];
        return;
    end
    
    try
        file_content = textscan(fid, '%s', 'Delimiter', '\n');
        lines = file_content{1};
        fclose(fid);
    catch
        fclose(fid);
        time_points = [];
        allX = [];
        allY = [];
        return;
    end
    
    % 查找ZONE
    zone_starts = find(contains(lines, 'ZONE T='));
    n_zones = length(zone_starts);
    n_time_steps = floor(n_zones / 2);
    
    % 预分配
    time_points = zeros(n_time_steps, 1);
    allX = cell(n_time_steps, 1);
    allY = cell(n_time_steps, 1);
    
    valid_steps = 0;
    
    % 处理每个时间步（第二个ZONE）
    for i = 1:n_time_steps
        flag_zone_idx = 2 * i;
        
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
        
        % 找到数据位置
        data_start = zone_starts(flag_zone_idx);
        while data_start <= length(lines) && ~contains(lines{data_start}, 'I=')
            data_start = data_start + 1;
        end
        data_start = data_start + 1;
        
        if flag_zone_idx < n_zones
            data_end = zone_starts(flag_zone_idx + 1) - 1;
        else
            data_end = length(lines);
        end
        
        % 读取数据
        [X_data, Y_data] = readZoneData(lines, data_start, data_end, expected_points);
        
        if length(X_data) >= expected_points && length(Y_data) >= expected_points
            valid_steps = valid_steps + 1;
            time_points(valid_steps) = current_time;
            allX{valid_steps} = X_data(1:expected_points);
            allY{valid_steps} = Y_data(1:expected_points);
        end
    end
    
    % 截取有效数据
    time_points = time_points(1:valid_steps);
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
    end
    
    X_data = X_data(:);
    Y_data = Y_data(:);
end

% 简化的曲率计算函数
function int_curv = calculateCurvature(allX, allY, time, mask)
    indices = find(mask);
    n_times = length(indices);
    int_curv = zeros(n_times, 1);
    
    for i = 1:n_times
        idx = indices(i);
        x = allX{idx};
        y = allY{idx};
        
        n_pts = length(x);
        if n_pts >= 3
            % 计算弧长
            ds = sqrt(diff(x).^2 + diff(y).^2);
            s = [0; cumsum(ds)];
            
            % 计算曲率
            kappa = zeros(n_pts, 1);
            
            % 内部点
            for j = 2:n_pts-1
                x_pts = x(j-1:j+1);
                y_pts = y(j-1:j+1);
                s_pts = s(j-1:j+1);
                
                if length(unique(s_pts)) == 3
                    px = polyfit(s_pts, x_pts, 2);
                    py = polyfit(s_pts, y_pts, 2);
                    
                    dx_ds = 2*px(1)*s(j) + px(2);
                    dy_ds = 2*py(1)*s(j) + py(2);
                    d2x_ds2 = 2*px(1);
                    d2y_ds2 = 2*py(1);
                    
                    numerator = dx_ds * d2y_ds2 - dy_ds * d2x_ds2;
                    denominator = (dx_ds^2 + dy_ds^2)^(3/2);
                    
                    if denominator > 1e-10
                        kappa(j) = numerator / denominator;
                    end
                end
            end
            
            % 边界点（简化处理）
            kappa(1) = kappa(2);
            kappa(end) = kappa(end-1);
            
            % 积分
            int_curv(i) = trapz(s, kappa);
        end
    end
end

fprintf('\n分析完成！\n');
