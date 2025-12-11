% 扭曲数绝对值随三元耦合强度变化的函数关系研究
% 计算不同扭曲数绝对值的占比随sigma_delta的变化函数
clc, clear, close all;

% 记录总开始时间
total_start_time = tic;

%% 参数设置
n = 83;                      % 振荡器数量
r = 2;                       % 耦合范围
sigma = 1;                   % 成对耦合强度（已通过时间缩放设为1）
sigma_delta_values = 0.2:0.6:7;  % 三元耦合强度取值范围
num_samples = 16000;          % 每个sigma_delta值的样本数量（可根据计算速度调整）

% 退火参数
T_anneal = 400;              % 退火总时间
dt_anneal = 0.1;             % 退火时间步长
anneal_steps = T_anneal/dt_anneal;
D_start = 5.0;               % 初始噪声强度（足够大以遍历整个空间）
D_end = 0.0;                 % 最终噪声强度

% 稳定参数
T_stable = 40;               % 退火后稳定时间
dt_stable = 0.1;             % 稳定阶段时间步长
stable_steps = T_stable/dt_stable;

% 总步数
total_steps = anneal_steps + stable_steps;

% 有序度阈值
order_threshold = 0.99;       % 有序度阈值，小于此值视为无效样本

%% 开始模拟
fprintf('参数设置:\n');
fprintf('振荡器数量: %d\n', n);
fprintf('耦合范围: %d\n', r);
fprintf('成对耦合强度: %.1f\n', sigma);
fprintf('三元耦合强度范围: %.1f:0.2:%.1f\n', sigma_delta_values(1), sigma_delta_values(end));
fprintf('样本数量: %d (每个sigma_delta值)\n', num_samples);
fprintf('初始噪声强度: %.1f\n', D_start);
fprintf('最终噪声强度: %.1f\n', D_end);
fprintf('退火时间: %.1f秒\n', T_anneal);
fprintf('稳定时间: %.1f秒\n', T_stable);
fprintf('有序度阈值: %.2f\n', order_threshold);
fprintf('使用随机初始状态\n\n');

% 初始化并行环境 - 使用全部32线程
if isempty(gcp('nocreate'))
    % 创建使用32个工作线程的本地并行池
    parpool('local'); % 利用全部16个大核心32个线程
end

% 创建噪声强度随时间的变化曲线（线性退火）
D_schedule = zeros(anneal_steps, 1);
for t = 1:anneal_steps
    % 线性降低噪声强度
    D_schedule(t) = D_start - (D_start - D_end) * (t-1)/(anneal_steps-1);
end

% 初始化结果数组
% 记录每个sigma_delta值下不同|q|的百分比和有效样本百分比
results = struct('sigma_delta', {}, 'abs_q_0', {}, 'abs_q_1', {}, 'abs_q_2', {}, ...
                 'abs_q_3', {}, 'abs_q_ge_4', {}, 'invalid', {}, 'valid_percent', {});

% 计算总共需要进行的sigma_delta模拟次数
total_delta_count = length(sigma_delta_values);
delta_progress = 0;

% 对每个sigma_delta值进行模拟
for delta_idx = 1:length(sigma_delta_values)
    % 显示整体进度
    delta_progress = delta_idx;
    total_progress = ((delta_idx-1) / total_delta_count * 100);
    fprintf('\n总体进度: %.1f%% [%d/%d] - ', total_progress, delta_idx, total_delta_count);
    
    % 当前三元耦合强度
    sigma_delta = sigma_delta_values(delta_idx);
    fprintf('模拟 σ_Δ = %.1f...\n', sigma_delta);
    
    % 记录此sigma_delta的开始时间
    sigma_delta_start_time = tic;
    
    % 初始化保存所有样本的扭曲数和有序度数据结构
    all_winding_numbers = zeros(num_samples, 1);
    all_orders = zeros(num_samples, 1);
    
    % 使用并行循环处理所有样本
    fprintf('开始并行量子退火模拟...\n');
    
    % 制作进度条变量
    progress_bar_length = 30;
    last_percent = 0;
    
    % 由于parfor不能直接显示进度，所以先进行模拟
    parfor sample_idx = 1:num_samples
        % 随机初始相位
        theta = 2*pi*rand(n, 1) - pi;  % 随机初始态 [-π, π]
        
        % 退火阶段 - 从高噪声到低噪声
        for t = 1:anneal_steps
            % 获取当前噪声强度
            D_current = D_schedule(t);
            
            % 计算相互作用力
            dtheta = compute_interactions(theta, n, r, sigma, sigma_delta);
            
            % 更新相位（带有当前噪声强度）
            noise = D_current * sqrt(2*dt_anneal) * randn(n, 1);
            theta = theta + dtheta * dt_anneal + noise;
            
            % 确保相位在[-π, π]范围内
            theta = mod(theta + pi, 2*pi) - pi;
        end
        
        % 稳定阶段 - 无噪声，让系统稳定
        for t = 1:stable_steps
            % 计算相互作用力
            dtheta = compute_interactions(theta, n, r, sigma, sigma_delta);
            
            % 更新相位（无噪声）
            theta = theta + dtheta * dt_stable;
            
            % 确保相位在[-π, π]范围内
            theta = mod(theta + pi, 2*pi) - pi;
        end
        
        % 计算最终扭曲数
        final_winding_number = calculate_winding_number(theta);
        
        % 计算最终状态的有序度
        final_order = calculate_order(theta);
        
        % 存储结果
        all_winding_numbers(sample_idx) = final_winding_number;
        all_orders(sample_idx) = final_order;
    end
    
    % 显示计算完成的进度条
    bar = repmat('#', 1, progress_bar_length);
    fprintf('\r计算完成: [%s] 100%%\n', bar);
    
    % 统计有序样本和不同扭曲数绝对值的比例
    valid_samples = all_orders >= order_threshold;
    valid_count = sum(valid_samples);
    valid_percent = (valid_count / num_samples) * 100;
    
    % 对有效样本计算不同|q|的比例
    if valid_count > 0
        valid_windings = all_winding_numbers(valid_samples);
        abs_valid_windings = abs(valid_windings);
        
        % 计算各|q|值的比例
        q_0_percent = sum(abs_valid_windings == 0) / valid_count * 100;
        q_1_percent = sum(abs_valid_windings == 1) / valid_count * 100;
        q_2_percent = sum(abs_valid_windings == 2) / valid_count * 100;
        q_3_percent = sum(abs_valid_windings == 3) / valid_count * 100;
        q_ge_4_percent = sum(abs_valid_windings >= 4) / valid_count * 100;
    else
        % 如果没有有效样本，所有比例都设为0
        q_0_percent = 0;
        q_1_percent = 0;
        q_2_percent = 0;
        q_3_percent = 0;
        q_ge_4_percent = 0;
    end
    
    invalid_percent = 100 - valid_percent;
    
    % 存储结果
    results(delta_idx).sigma_delta = sigma_delta;
    results(delta_idx).abs_q_0 = q_0_percent;
    results(delta_idx).abs_q_1 = q_1_percent;
    results(delta_idx).abs_q_2 = q_2_percent;
    results(delta_idx).abs_q_3 = q_3_percent;
    results(delta_idx).abs_q_ge_4 = q_ge_4_percent;
    results(delta_idx).invalid = invalid_percent;
    results(delta_idx).valid_percent = valid_percent;
    
    % 计算并显示此sigma_delta的运行时间
    sigma_delta_time = toc(sigma_delta_start_time);
    fprintf('σ_Δ = %.1f 计算完成，用时 %.2f 秒 (%.2f 分钟)\n', sigma_delta, sigma_delta_time, sigma_delta_time/60);
    fprintf('有效样本比例: %.2f%%\n', valid_percent);
    
    % 显示总体进度
    total_progress = (delta_idx / total_delta_count * 100);
    fprintf('总体进度: %.1f%% [%d/%d]\n', total_progress, delta_idx, total_delta_count);
end

% 计算并显示总运行时间
total_time = toc(total_start_time);
fprintf('\n所有模拟完成，总运行时间: %.2f秒 (%.2f分钟)\n', total_time, total_time/60);

%% 准备绘图数据
sigma_delta_range = [results.sigma_delta];
abs_q_0_data = [results.abs_q_0];
abs_q_1_data = [results.abs_q_1];
abs_q_2_data = [results.abs_q_2];
abs_q_3_data = [results.abs_q_3];
abs_q_ge_4_data = [results.abs_q_ge_4];
invalid_data = [results.invalid];

%% 设置图表通用样式
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);
set(0, 'DefaultLineLineWidth', 2.0);

%% 绘制|q|随sigma_delta变化的函数图
figure('Position', [100, 100, 1000, 700], 'Color', 'white');

% 使用科学可视化专业配色方案
scientific_colors = [
    [31, 119, 180]/255;    % 蓝色 (|q| = 0)
    [44, 160, 44]/255;     % 绿色 (|q| = 1)
    [148, 103, 189]/255;   % 紫色 (|q| = 2)
    [227, 119, 194]/255;   % 粉色 (|q| = 3)
    [214, 39, 40]/255;     % 红色 (|q| >= 4)
    [127, 127, 127]/255    % 灰色 (无效样本)
];

% 绘制数据点
hold on;
h1 = scatter(sigma_delta_range, abs_q_0_data, 80, scientific_colors(1,:), 'filled', 'DisplayName', '|q| = 0');
h2 = scatter(sigma_delta_range, abs_q_1_data, 80, scientific_colors(2,:), 'filled', 'DisplayName', '|q| = 1');
h3 = scatter(sigma_delta_range, abs_q_2_data, 80, scientific_colors(3,:), 'filled', 'DisplayName', '|q| = 2');
h4 = scatter(sigma_delta_range, abs_q_3_data, 80, scientific_colors(4,:), 'filled', 'DisplayName', '|q| = 3');
h5 = scatter(sigma_delta_range, abs_q_ge_4_data, 80, scientific_colors(5,:), 'filled', 'DisplayName', '|q| \geq 4');
h6 = scatter(sigma_delta_range, invalid_data, 80, scientific_colors(6,:), 'filled', 'DisplayName', 'Invalid samples');

% 拟合曲线部分 - 使用更平滑的样条插值
if length(sigma_delta_range) >= 5
    finer_sigma_delta = linspace(sigma_delta_values(1), sigma_delta_values(end), 500); % 增加采样点数量，使曲线更平滑
    
    % 使用形状保持的插值方法(PCHIP)或更平滑的样条插值
    % 根据您的图看起来更适合使用平滑样条而不是PCHIP
    abs_q_0_fit = spline(sigma_delta_range, abs_q_0_data, finer_sigma_delta);
    abs_q_1_fit = spline(sigma_delta_range, abs_q_1_data, finer_sigma_delta);
    abs_q_2_fit = spline(sigma_delta_range, abs_q_2_data, finer_sigma_delta);
    abs_q_3_fit = spline(sigma_delta_range, abs_q_3_data, finer_sigma_delta);
    abs_q_ge_4_fit = spline(sigma_delta_range, abs_q_ge_4_data, finer_sigma_delta);
    invalid_fit = spline(sigma_delta_range, invalid_data, finer_sigma_delta);
    
    % 使用较宽的线条增强视觉效果
    plot(finer_sigma_delta, abs_q_0_fit, 'Color', scientific_colors(1,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot(finer_sigma_delta, abs_q_1_fit, 'Color', scientific_colors(2,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot(finer_sigma_delta, abs_q_2_fit, 'Color', scientific_colors(3,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot(finer_sigma_delta, abs_q_3_fit, 'Color', scientific_colors(4,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot(finer_sigma_delta, abs_q_ge_4_fit, 'Color', scientific_colors(5,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot(finer_sigma_delta, invalid_fit, 'Color', scientific_colors(6,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
else
    % 直接连线部分
    plot(sigma_delta_range, abs_q_0_data, 'Color', scientific_colors(1,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot(sigma_delta_range, abs_q_1_data, 'Color', scientific_colors(2,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot(sigma_delta_range, abs_q_2_data, 'Color', scientific_colors(3,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot(sigma_delta_range, abs_q_3_data, 'Color', scientific_colors(4,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot(sigma_delta_range, abs_q_ge_4_data, 'Color', scientific_colors(5,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
    plot(sigma_delta_range, invalid_data, 'Color', scientific_colors(6,:), 'LineWidth', 2.5, 'HandleVisibility', 'off');
end
hold off;

% 添加网格和图例
grid on;
legend('Location', 'eastoutside', 'FontSize', 12);
xlabel('Three-Body Coupling Strength (σ_{\Delta})', 'FontSize', 16);
ylabel('Percentage (%)', 'FontSize', 16);
title('Winding Number Distribution vs. Three-Body Coupling Strength', 'FontSize', 20, 'FontWeight', 'bold');
subtitle(sprintf('P_{order} threshold = %.2f, %d samples per point', order_threshold, num_samples), 'FontSize', 14);

% 设置坐标轴范围和样式
xlim([sigma_delta_values(1) - 0.2, sigma_delta_values(end) + 0.2]);
ylim([0, 50]);
set(gca, 'Box', 'on', 'LineWidth', 1.5);
% 设置更精细的网格间隔
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XGrid', 'on', 'YGrid', 'on');
set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', ':');
set(gca, 'GridAlpha', 0.15, 'MinorGridAlpha', 0.1);


% 保存结果数据
save('winding_number_vs_sigma_delta3.mat', 'results', 'sigma_delta_range', 'abs_q_0_data', 'abs_q_1_data', ...
     'abs_q_2_data', 'abs_q_3_data', 'abs_q_ge_4_data', 'invalid_data');

% 输出表格格式的结果
fprintf('\n结果摘要（百分比 %%）:\n');
fprintf('σ_Δ\t|q|=0\t|q|=1\t|q|=2\t|q|=3\t|q|≥4\t无效样本\t有效样本比例\n');
fprintf('------------------------------------------------------------------------\n');
for i = 1:length(results)
    fprintf('%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t%.1f\t\t%.1f\n', ...
        results(i).sigma_delta, results(i).abs_q_0, results(i).abs_q_1, ...
        results(i).abs_q_2, results(i).abs_q_3, results(i).abs_q_ge_4, ...
        results(i).invalid, results(i).valid_percent);
end

%% 辅助函数1：计算相互作用
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

%% 辅助函数2：计算扭曲数
function winding_number = calculate_winding_number(theta)
    n = length(theta);
    phase_diffs = zeros(n, 1);
    
    for i = 1:n-1
        phase_diffs(i) = mod(theta(i+1) - theta(i) + pi, 2*pi) - pi;
    end
    phase_diffs(n) = mod(theta(1) - theta(n) + pi, 2*pi) - pi;
    
    % 计算扭曲数
    winding_number = round(sum(phase_diffs) / (2*pi));
end

%% 辅助函数3：计算有序度 - 局部有序振荡器比例
function p_order = calculate_order(theta)
    % 使用系统参数
    n = length(theta);
    r = 2;  % 耦合范围，与主程序保持一致
    
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