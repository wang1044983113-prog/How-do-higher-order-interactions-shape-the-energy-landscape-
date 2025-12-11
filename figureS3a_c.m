
clc, clear, close all;

% 记录总体开始时间
total_tic = tic;

% 检查并行计算工具箱是否可用，尝试创建并行池
try
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool; % 创建默认大小的并行池
    end
    fprintf('并行计算已启用，使用 %d 个工作线程\n', poolobj.NumWorkers);
catch
    warning('无法启用并行计算，请确保已安装Parallel Computing Toolbox');
end

% 参数设置
n = 83;          % 振荡器数量
r = 2;           % 耦合范围
sigma = 1;       % 成对耦合强度

% 多个三元耦合强度比较
sigma_delta_values = [2.5,3,3.5];
num_sigma_delta = length(sigma_delta_values);

prep_time = 100;     % 无噪声状态准备时间（确定初始状态）
max_sim_time = 2000; % 最大模拟时间（计算离出时间）
dt = 0.1;           % 时间步长
prep_steps = prep_time/dt;
max_steps = max_sim_time/dt;

% 不同的去噪声时间参数
no_noise_time_values = [2,5,10];    % 去掉噪声后计算的时间
num_no_noise_times = length(no_noise_time_values);
exit_threshold = 0.999;  % 固定离出阈值

% 检测间隔设置 - 优化检测频率
detect_interval = 1;  % 每1秒检测一次状态，而不是每一步
detect_steps = round(detect_interval/dt);  % 对应的步数间隔

% 初始non-twisted阈值（仍然使用固定值来生成初始样本）
initial_nontwisted_threshold = 0.8;  % P_order < 0.8 为初始non-twisted状态

% 设置噪声强度范围 - 使用均匀分布的1/D²值重新计算D_values
inv_D_squared_min = 1/(0.4^2);  % 对应于D=0.4的1/D²
inv_D_squared_max = 1/(0.25^2);  % 对应于D=0.25的1/D²
num_D_values = 10;  % 设置点数

% 创建均匀分布的1/D²值
inv_D_squared_values = linspace(inv_D_squared_min, inv_D_squared_max, num_D_values);

% 计算对应的D值
D_values = 1 ./ sqrt(inv_D_squared_values);

% 样本数量
num_samples = 1600;  % 每个状态类型的样本数量

% 存储不同no_noise_time和sigma_delta值的结果
mean_exit_time_all = zeros(num_no_noise_times, num_sigma_delta, length(D_values));  % 平均MFPT
std_exit_time_all = zeros(num_no_noise_times, num_sigma_delta, length(D_values));  % 标准差
valid_samples_all = zeros(num_no_noise_times, num_sigma_delta, length(D_values));  % 有效样本数
all_exit_times_data = cell(num_no_noise_times, num_sigma_delta, length(D_values));  % 存储所有原始有效离出时间数据

% 先生成初始状态样本（对所有情况共用）
fprintf('第一阶段：生成初始Non-twisted状态样本...\n');
prep_tic = tic;

% 生成non-twisted的初始状态样本
nontwisted_samples = cell(num_samples, 1);

% 使用并行计算生成样本
result_cell = cell(num_samples * 6, 1);
parfor sample_idx = 1:num_samples * 6  % 生成更多的样本，以确保获得足够的non-twisted状态
    % 随机初始相位
    initial_phases = 2*pi*rand(n, 1) - pi;
    
    % 无噪声系统演化，确定初始状态
    [state_type, theta_final] = determine_initial_state(initial_phases, n, r, sigma, sigma_delta_values(1), prep_steps, dt, initial_nontwisted_threshold);
    
    % 返回状态和最终相位
    result = struct('state_type', state_type, 'theta_final', theta_final);
    result.sample_idx = sample_idx;
    
    % 标记类别
    if strcmp(state_type, 'twisted')
        result.category = 'twisted';
    else
        result.category = 'nontwisted';
    end
    
    % 返回结果
    result_cell{sample_idx} = result;
end

% 在并行计算完成后，分配样本
valid_results = [];
for i = 1:length(result_cell)
    if ~isempty(result_cell{i})
        valid_results = [valid_results; result_cell{i}];
    end
end

% 提取non-twisted状态的样本
nontwisted_results = valid_results(strcmp({valid_results.category}, 'nontwisted'));

% 确保我们有足够的non-twisted样本
if length(nontwisted_results) < num_samples
    error('没有生成足够的Non-twisted样本。获得: %d, 需要: %d', length(nontwisted_results), num_samples);
end

% 只取所需数量的样本
nontwisted_samples = nontwisted_results(1:num_samples);

prep_time_elapsed = toc(prep_tic);
fprintf('样本生成完成，用时 %.2f 秒。生成 %d 个non-twisted样本。\n', prep_time_elapsed, num_samples);

% 对每个no_noise_time值和sigma_delta值进行模拟
for nt_idx = 1:num_no_noise_times
    no_noise_time = no_noise_time_values(nt_idx);
    no_noise_steps = no_noise_time/dt;
    
    fprintf('\n\n===== 开始计算 no_noise_time = %d 秒的情况 (%d/%d) =====\n\n', ...
            no_noise_time, nt_idx, num_no_noise_times);
    
    for sd_idx = 1:num_sigma_delta
        sigma_delta = sigma_delta_values(sd_idx);
        fprintf('\n--- sigma_delta = %.1f (%d/%d) ---\n', ...
                sigma_delta, sd_idx, num_sigma_delta);
        
        % 存储平均首次离出时间数据
        mean_exit_time = zeros(size(D_values));
        std_exit_time = zeros(size(D_values));
        total_samples = num_samples;  % 已知总样本数
        valid_samples = zeros(size(D_values));
        
        fprintf('计算Non-twisted到Twisted的离出时间 (噪声强度 %.2f-%.2f):\n', min(D_values), max(D_values));
        fprintf('离出条件：每个时间点去掉噪声，计算 %d 秒后的P_order，大于 %.2f 算离出\n', no_noise_time, exit_threshold);
        
        % 对每个D值进行处理
        for d_idx = 1:length(D_values)
            % 记录当前D值的计算开始时间
            d_tic = tic;
            
            D = D_values(d_idx);
            fprintf('  处理噪声强度 D = %.4f (%d/%d)\n', D, d_idx, length(D_values));
            
            % 计算离出时间
            all_exit_times = zeros(num_samples, 1);
            
            parfor sample_idx = 1:num_samples
                % 获取预先计算的non-twisted样本
                sample = nontwisted_samples(sample_idx);
                theta_final = sample.theta_final;
                
                % 计算从non-twisted状态到twisted状态的离出时间
                exit_time = calculate_exit_time_new_condition(theta_final, n, r, sigma, sigma_delta, D, max_steps, dt, no_noise_steps, exit_threshold, detect_steps);
                
                % 存储结果
                all_exit_times(sample_idx) = exit_time;
            end
            
            % 处理结果 - 只考虑有效离出的样本
            valid_indices = all_exit_times < max_sim_time;
            valid_exit_times = all_exit_times(valid_indices);
            valid_samples(d_idx) = sum(valid_indices);
            
            % 保存原始有效离出时间数据用于后续误差棒分析
            all_exit_times_data{nt_idx, sd_idx, d_idx} = valid_exit_times;
            
            % 计算平均离出时间和标准差 - 只使用有效样本（方法1：先求平均）
            if valid_samples(d_idx) > 0
                mean_exit_time(d_idx) = mean(valid_exit_times);
                std_exit_time(d_idx) = std(valid_exit_times);
            else
                % 如果没有有效样本，设置为NaN
                mean_exit_time(d_idx) = NaN;
                std_exit_time(d_idx) = NaN;
            end
            
            % 记录当前D值的计算时间
            compute_time_per_D = toc(d_tic);
            
            % 输出当前进度和统计信息
            if valid_samples(d_idx) > 0
                fprintf('    %d/%d (有效样本, %.1f%%), MFPT = %.2f, 用时 %.1fs\n', ...
                    valid_samples(d_idx), total_samples, ...
                    100*valid_samples(d_idx)/total_samples, mean_exit_time(d_idx), compute_time_per_D);
            else
                fprintf('    %d/%d (有效样本, %.1f%%), MFPT = NaN, 用时 %.1fs\n', ...
                    valid_samples(d_idx), total_samples, ...
                    100*valid_samples(d_idx)/total_samples, compute_time_per_D);
            end
        end
        
        % 存储当前no_noise_time和sigma_delta的结果
        mean_exit_time_all(nt_idx, sd_idx, :) = mean_exit_time;
        std_exit_time_all(nt_idx, sd_idx, :) = std_exit_time;
        valid_samples_all(nt_idx, sd_idx, :) = valid_samples;
    end
end

% 记录总计算时间
total_compute_time = toc(total_tic);
fprintf('\n总计算时间: %.2f 秒 (%.2f 分钟, %.2f 小时)\n', ...
    total_compute_time, total_compute_time/60, total_compute_time/3600);

%% 绘制三个子图
fprintf('\n第三阶段：绘制三个子图的ln(MFPT)与1/D²关系图\n');

% 创建新的横坐标：1/D²
x_new = 1 ./ (D_values.^2);

% 计算主数据点：使用方法1（先求平均再取对数）
y_new_all = zeros(num_no_noise_times, num_sigma_delta, length(D_values));
for nt_idx = 1:num_no_noise_times
    for sd_idx = 1:num_sigma_delta
        valid_indices = ~isnan(squeeze(mean_exit_time_all(nt_idx, sd_idx, :)));
        temp_data = squeeze(mean_exit_time_all(nt_idx, sd_idx, :));
        y_new_all(nt_idx, sd_idx, valid_indices) = log(temp_data(valid_indices));  % 先平均再取对数
    end
end

% 计算误差棒：使用方法2（先取对数再求平均）
y_err_all = zeros(num_no_noise_times, num_sigma_delta, length(D_values));  % 对数尺度标准差

for nt_idx = 1:num_no_noise_times
    for sd_idx = 1:num_sigma_delta
        for i = 1:length(D_values)
            if valid_samples_all(nt_idx, sd_idx, i) > 0 && ~isempty(all_exit_times_data{nt_idx, sd_idx, i})
                % 获取原始有效离出时间数据
                valid_exit_times = all_exit_times_data{nt_idx, sd_idx, i};
                
                % 先对每个数据点取对数
                log_exit_times = log(valid_exit_times);
                
                % 然后计算对数数据的标准差
                y_err_all(nt_idx, sd_idx, i) = std(log_exit_times);
            else
                y_err_all(nt_idx, sd_idx, i) = NaN;
            end
        end
    end
end

% 创建颜色方案 - 为每个sigma_delta值分配不同颜色
colors = [
    0.8500, 0.3250, 0.0980;  % 红橙色
    0.0, 0.4470, 0.7410;     % 蓝色
    0.4660, 0.6740, 0.1880;  % 绿色
];

% 创建适中的背景颜色（用于柱状图）
light_colors = colors + 0.25 * (1 - colors);  % 向白色方向混合

% 设置线型和标记类型
marker_styles = {'o', 's', '^'};
marker_sizes = [12, 12, 12];  % 标记大小

% 计算拟合参数
p_values = zeros(num_no_noise_times, num_sigma_delta, 2);  % 存储拟合参数
R_squared_values = zeros(num_no_noise_times, num_sigma_delta);  % 存储R²值

for nt_idx = 1:num_no_noise_times
    for sd_idx = 1:num_sigma_delta
        % 提取当前条件下的数据
        y_temp = squeeze(y_new_all(nt_idx, sd_idx, :));
        valid_indices = ~isnan(y_temp);
        x_valid = x_new(valid_indices);
        y_valid = y_temp(valid_indices);
        
        if length(x_valid) >= 2 && length(y_valid) >= 2
            % 确保x_valid和y_valid都是列向量
            x_valid = x_valid(:);
            y_valid = y_valid(:);
            
            % 线性拟合: ln(MFPT) = a + b*(1/D²)
            p = polyfit(x_valid, y_valid, 1);
            p_values(nt_idx, sd_idx, :) = p;
            
            % 计算R²值
            y_pred = polyval(p, x_valid);
            SSE = sum((y_valid - y_pred).^2);
            SST = sum((y_valid - mean(y_valid)).^2);
            if SST > 0
                R_squared = 1 - SSE/SST;
            else
                R_squared = NaN;
            end
            R_squared_values(nt_idx, sd_idx) = R_squared;
        else
            % 如果数据不足，填入NaN
            p_values(nt_idx, sd_idx, :) = [NaN, NaN];
            R_squared_values(nt_idx, sd_idx) = NaN;
        end
    end
end

% 创建包含三个子图的大图
fig1 = figure('Position', [50, 50, 1800, 600], 'Color', 'white');
set(gcf, 'DefaultAxesFontName', 'Arial');

subplot_titles = {
    sprintf('no\\_noise\\_time = %d s', no_noise_time_values(1)),
    sprintf('no\\_noise\\_time = %d s', no_noise_time_values(2)),
    sprintf('no\\_noise\\_time = %d s', no_noise_time_values(3))
};

% 绘制三个子图
for nt_idx = 1:num_no_noise_times
    subplot(1, 3, nt_idx);
    hold on;
    
    % 第一步：创建左坐标轴并绘制主要数据（前景层）
    yyaxis left;
    hold on;
    
    h_data = zeros(num_sigma_delta, 1);
    h_plot = zeros(num_sigma_delta, 1);
    
    % 先绘制拟合线
    for sd_idx = 1:num_sigma_delta
        % 提取当前条件下的数据
        y_temp = squeeze(y_new_all(nt_idx, sd_idx, :));
        valid_indices = ~isnan(y_temp);
        x_valid = x_new(valid_indices);
        y_valid = y_temp(valid_indices);
        
        if length(x_valid) >= 2 && ~any(isnan(squeeze(p_values(nt_idx, sd_idx, :))))
            x_fit = linspace(6, 16.3, 100);  % 固定横坐标范围
            y_fit = polyval(squeeze(p_values(nt_idx, sd_idx, :)), x_fit);
            
            h_plot(sd_idx) = plot(x_fit, y_fit, '-', 'LineWidth', 2.5, 'Color', colors(sd_idx, :));
        end
    end
    
    % 再绘制数据点（最前景层）
    for sd_idx = 1:num_sigma_delta
        % 提取当前条件下的数据
        y_temp = squeeze(y_new_all(nt_idx, sd_idx, :));
        valid_indices = ~isnan(y_temp);
        x_valid = x_new(valid_indices);
        y_valid = y_temp(valid_indices);
        
        if ~isempty(x_valid) && ~isempty(y_valid)
            h_data(sd_idx) = plot(x_valid, y_valid, marker_styles{sd_idx}, ...
                'MarkerSize', marker_sizes(sd_idx), ...
                'MarkerFaceColor', colors(sd_idx, :), ...
                'MarkerEdgeColor', 'none', ...
                'Color', colors(sd_idx, :), ...
                'LineWidth', 2.0, ...
                'LineStyle', 'none');
        end
    end
    
    % 设置左坐标轴
    if nt_idx == 1
        ylabel('ln(\tau_e)', 'FontSize', 20, 'FontWeight', 'bold');
    end
    set(gca, 'YColor', 'k', 'FontSize', 14, 'FontWeight', 'bold');
    
    % 第二步：切换到右坐标轴绘制误差柱状图（背景层）
    yyaxis right;
    hold on;
    
    % 计算柱状图位置 - 调整柱子宽度
    bar_width = 0.18;  % 调整柱子宽度适应横坐标范围
    group_width = bar_width * num_sigma_delta;
    
    % 绘制所有柱状图（背景层，会自动在主图下方）
    for i = 1:length(x_new)  % 所有x位置
        start_pos = x_new(i) - group_width/2 + bar_width/2;
        for sd_idx = 1:num_sigma_delta  % 每个位置的柱子
            if ~isnan(y_err_all(nt_idx, sd_idx, i))
                bar_pos = start_pos + (sd_idx-1) * bar_width;
                bar(bar_pos, y_err_all(nt_idx, sd_idx, i), bar_width, ...
                    'FaceColor', light_colors(sd_idx, :), ...
                    'EdgeColor', 'none', ...
                    'FaceAlpha', 0.25);  % 透明的背景
            end
        end
    end
    
    % 设置右坐标轴
    if nt_idx == 3
        ylabel('Standard Deviation of ln(\tau_e)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', [0.0, 0.3, 0.7]);
    end
    set(gca, 'YColor', [0.0, 0.3, 0.7], 'FontSize', 14, 'FontWeight', 'bold');
    
    % 第三步：重新激活左坐标轴进行最终设置
    yyaxis left;
    
    % 设置x轴
    xlabel('1/D^2', 'FontSize', 20, 'FontWeight', 'bold');
    xlim([5.8, 16.5]);

    % 设置图表标题
title(sprintf('Remove noise for %d seconds', no_noise_time_values(nt_idx)), 'FontSize', 20, 'FontWeight', 'bold');
    
    % 设置固定的y轴范围
    ylim([3, 6]);  % 左侧y轴范围设置为3-6
    
    % 第四步：设置右坐标轴的y轴范围
    yyaxis right;
    ylim([0, 6]);  % 右侧y轴范围设置为0-6
    
    % 第五步：重新激活左坐标轴完成其他设置
    yyaxis left;
    
    % 设置图表外观
    box on;
    set(gca, 'FontSize', 14, 'FontWeight', 'bold', ...
        'LineWidth', 1.0, ...
        'TickDir', 'out');
    


    % 为每个子图添加图例，显示拟合的斜率值（k值）
    legend_handles = [];
    legend_str = {};
    legend_count = 0;
    
    for sd_idx = 1:num_sigma_delta
        if ishandle(h_data(sd_idx)) && h_data(sd_idx) ~= 0
            legend_count = legend_count + 1;
            legend_handles(legend_count) = h_data(sd_idx);
            % 所有子图都显示详细信息，包括拟合的斜率值
            if ~isnan(p_values(nt_idx, sd_idx, 1))
                legend_str{legend_count} = sprintf('\\sigma_{\\Delta} = %.1f (k ≈ %.3f)', ...
                    sigma_delta_values(sd_idx), p_values(nt_idx, sd_idx, 1));
            else
                legend_str{legend_count} = sprintf('\\sigma_{\\Delta} = %.1f (no fit)', ...
                    sigma_delta_values(sd_idx));
            end
        end
    end
    
    % 只有当有有效句柄时才添加图例
    if ~isempty(legend_handles)
        legend(legend_handles, legend_str, 'Location', 'northwest', 'FontSize', 10, 'FontWeight', 'bold');
    end
end

% 添加整体标题
sgtitle('ln(\tau_e) vs. 1/D^2 for Different no\_noise\_time Values', 'FontSize', 18, 'FontWeight', 'bold');

% % 保存所有结果数据
% save('mean_exit_time_results_three_subplots.mat', ...
%     'D_values', 'sigma_delta_values', 'no_noise_time_values', 'mean_exit_time_all', 'std_exit_time_all', ...
%     'valid_samples_all', 'all_exit_times_data', 'p_values', 'R_squared_values', 'y_err_all', ...
%     'exit_threshold', ...
%     'n', 'r', 'sigma', 'num_samples', 'prep_time', 'max_sim_time', ...
%     'initial_nontwisted_threshold', 'total_compute_time');
% 
% % 保存图片
% saveas(fig1, 'ln_tau_vs_1_D2_three_subplots.png');

% 显示拟合结果信息
fprintf('\n不同no_noise_time和sigma_delta值的拟合结果对比：\n');
fprintf('=======================================================\n');
for nt_idx = 1:num_no_noise_times
    fprintf('no_noise_time = %d s:\n', no_noise_time_values(nt_idx));
    fprintf('-------------------------------------------------------\n');
    for sd_idx = 1:num_sigma_delta
        fprintf('  sigma_delta = %.1f:\n', sigma_delta_values(sd_idx));
        if ~isnan(p_values(nt_idx, sd_idx, 1))
            fprintf('    ln(MFPT) = %.4f + %.4f·(1/D²)\n', p_values(nt_idx, sd_idx, 2), p_values(nt_idx, sd_idx, 1));
            fprintf('    R² = %.4f\n', R_squared_values(nt_idx, sd_idx));
            fprintf('    能量势垒高度ΔV ≈ %.4f\n', p_values(nt_idx, sd_idx, 1));
        else
            fprintf('    无法进行拟合（数据不足）\n');
        end
    end
    fprintf('-------------------------------------------------------\n');
end

fprintf('\n混合方法说明:\n');
fprintf('- 主数据点: 先求平均离出时间，再取对数 ln(mean(τ))\n');
fprintf('- 误差棒: 先对每个数据点取对数，再求标准差 std(ln(τ))\n');
fprintf('已生成三个子图并保存\n');

%% ============= 辅助函数 =============

% 辅助函数1: 确定初始状态（无噪声系统演化）
function [state_type, final_theta] = determine_initial_state(initial_phases, n, r, sigma, sigma_delta, steps, dt, p_order_threshold)
    theta = initial_phases;
    
    % 无噪声系统演化
    for t = 1:steps
        dtheta = compute_interactions(theta, n, r, sigma, sigma_delta);
        theta = theta + dtheta * dt;
        theta = mod(theta + pi, 2*pi) - pi;
    end
    
    % 计算P_order
    p_order = calculate_p_order(theta, n, r);
    
    % 确定状态类型
    if p_order > p_order_threshold
        state_type = 'twisted';
    else
        state_type = 'nontwisted';
    end
    
    final_theta = theta;
end
% 辅助函数2: 计算首次离出时间（修改后的离出条件）
function exit_time = calculate_exit_time_new_condition(initial_theta, n, r, sigma, sigma_delta, D, max_steps, dt, no_noise_steps, exit_threshold, detect_steps)
    theta = initial_theta;
    exit_time = max_steps * dt;  % 默认为最大模拟时间
    
    % 添加噪声的系统演化
    t = 1;
    while t <= max_steps
        % 向前模拟detect_steps步或直到最大步数
        for i = 1:min(detect_steps, max_steps-t+1)
            dtheta = compute_interactions(theta, n, r, sigma, sigma_delta);
            noise = D * sqrt(2 * dt) * randn(n, 1);
            theta = theta + dtheta * dt + noise;
            theta = mod(theta + pi, 2*pi) - pi;
        end
        
        t = t + detect_steps;
        
        % 新离出条件：复制当前状态，模拟无噪声演化指定时间
        theta_no_noise = theta;  % 复制当前状态
        
        % 无噪声演化指定时间，同时记录p_order
        no_noise_start_step = round(1/dt);  % 从1秒开始记录
        p_order_values = [];  % 存储p_order值
        
        for i = 1:no_noise_steps
            dtheta = compute_interactions(theta_no_noise, n, r, sigma, sigma_delta);
            theta_no_noise = theta_no_noise + dtheta * dt;  % 无噪声
            theta_no_noise = mod(theta_no_noise + pi, 2*pi) - pi;
            
            % 从1秒开始记录p_order
            if i >= no_noise_start_step
                p_order_current = calculate_p_order(theta_no_noise, n, r);
                p_order_values = [p_order_values; p_order_current];
            end
        end
        
        % 计算平均p_order（从1秒到设置时间结束）
        if ~isempty(p_order_values)
            mean_p_order = mean(p_order_values);
        else
            mean_p_order = calculate_p_order(theta_no_noise, n, r);  % 如果时间太短，只用最后值
        end
        
        % 检查是否满足离出条件
        if mean_p_order > exit_threshold
            exit_time = (t - detect_steps) * dt;  % 记录离出时间
            break;
        end
    end
    
    return;
end

% 辅助函数3: 计算相互作用
function dtheta = compute_interactions(theta, n, r, sigma, sigma_delta)
    dtheta = zeros(n, 1);
    
    % 计算所有振荡器的相位变化
    for i = 1:n
        % 计算成对相互作用项
        pair_interaction = 0;
        for j = i-r:i+r
            % 处理周期性边界条件
            j_idx = mod(j-1, n) + 1;
            if j_idx ~= i
                pair_interaction = pair_interaction + sin(theta(j_idx) - theta(i));
            end
        end
        pair_interaction = sigma * pair_interaction / (2*r);
        
        % 计算三元相互作用项
        triad_interaction = 0;
        count = 0;
        for j = i-r:i+r
            for k = i-r:i+r
                % 处理周期性边界条件
                j_idx = mod(j-1, n) + 1;
                k_idx = mod(k-1, n) + 1;
                % 确保三角形涉及三个不同的节点
                if j_idx ~= i && k_idx ~= i && j_idx ~= k_idx
                    triad_interaction = triad_interaction + sin(theta(j_idx) + theta(k_idx) - 2*theta(i));
                    count = count + 1;
                end
            end
        end
        if count > 0
            triad_interaction = sigma_delta * triad_interaction / (2*r*(2*r-1));
        end
        
        % 合并两种相互作用
        dtheta(i) = pair_interaction + triad_interaction;
    end
end

% 辅助函数4: 计算P_order
function p_order = calculate_p_order(theta, n, r)
    % 计算局部序参量
    local_order = zeros(n, 1);
    for j = 1:n
        % 计算局部平均
        local_sum = 0;
        for k = j-r:j+r
            % 处理周期性边界条件
            k_idx = mod(k-1, n) + 1;
            local_sum = local_sum + exp(1i*theta(k_idx));
        end
        local_order(j) = abs(local_sum/(2*r+1));
    end
    
    % 计算有序振荡器比例
    disordered_count = sum(local_order < 0.85);
    p_order = 1 - disordered_count/n;
end