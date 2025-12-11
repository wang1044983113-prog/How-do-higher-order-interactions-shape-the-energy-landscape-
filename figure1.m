% 计算不同相位状态占比随sigma_delta的变化函数
clc, clear, close all;

% 记录总开始时间
total_start_time = tic;

%% 参数设置
n = 83;                      % 振荡器数量
r = 2;                       % 耦合范围
sigma = 1;                   % 成对耦合强度（已通过时间缩放设为1）
sigma_delta_values = 0:0.6:6;  % 三元耦合强度取值范围
num_samples = 160000;            % 每个sigma_delta值的样本数量

% 无噪声模拟参数
T_sim = 200;                 % 总模拟时间（200秒）
dt_sim = 0.1;                % 时间步长
sim_steps = T_sim/dt_sim;    % 总步数

% 有序度阈值
order_threshold = 0.99;      % 有序度阈值，大于等于此值视为twisted态

%% 开始模拟
fprintf('参数设置:\n');
fprintf('振荡器数量: %d\n', n);
fprintf('耦合范围: %d\n', r);
fprintf('成对耦合强度: %.1f\n', sigma);
fprintf('三元耦合强度范围: %.1f:0.6:%.1f\n', sigma_delta_values(1), sigma_delta_values(end));
fprintf('样本数量: %d (每个sigma_delta值)\n', num_samples);
fprintf('模拟时间: %.1f秒\n', T_sim);
fprintf('时间步长: %.1f\n', dt_sim);
fprintf('有序度阈值: %.2f\n', order_threshold);
fprintf('无噪声模拟\n\n');

% 初始化并行环境
if isempty(gcp('nocreate'))
    parpool('local'); % 利用全部可用线程
end

% 初始化结果数组
results = struct('sigma_delta', {}, 'twisted', {}, 'nontwisted', {});

% 计算总共需要进行的sigma_delta模拟次数
total_delta_count = length(sigma_delta_values);

% 对每个sigma_delta值进行模拟
for delta_idx = 1:length(sigma_delta_values)
    % 当前三元耦合强度
    sigma_delta = sigma_delta_values(delta_idx);
    fprintf('模拟 σ_Δ = %.1f [%d/%d]...\n', sigma_delta, delta_idx, total_delta_count);
    
    % 记录此sigma_delta的开始时间
    sigma_delta_start_time = tic;
    
    % 初始化保存所有样本的有序度数据结构
    all_orders = zeros(num_samples, 1);
    
    % 使用并行循环处理所有样本
    fprintf('开始并行无噪声模拟 (σ_Δ = %.1f)...\n', sigma_delta);
    
    % 显示进度的设置
    batch_size = ceil(num_samples / 100); % 将进度分成大约100个批次
    progress_bar_length = 50;
    
    for batch = 1:ceil(num_samples/batch_size)
        % 计算当前批次的样本范围
        start_idx = (batch-1) * batch_size + 1;
        end_idx = min(batch * batch_size, num_samples);
        current_batch_size = end_idx - start_idx + 1;
        
        % 创建当前批次的临时结果数组
        batch_orders = zeros(current_batch_size, 1);
        
        % 使用parfor处理当前批次的样本
        parfor i = 1:current_batch_size
            % 随机初始相位
            theta = 2*pi*rand(n, 1) - pi;  % 随机初始态 [-π, π]
            
            % 无噪声模拟
            for t = 1:sim_steps
                % 计算相互作用力
                dtheta = compute_interactions(theta, n, r, sigma, sigma_delta);
                
                % 更新相位（无噪声）
                theta = theta + dtheta * dt_sim;
                
                % 确保相位在[-π, π]范围内
                theta = mod(theta + pi, 2*pi) - pi;
            end
            
            % 计算最终状态的有序度
            final_order = calculate_order(theta);
            
            % 存储结果到临时数组
            batch_orders(i) = final_order;
        end
        
        % 将批次结果复制到主结果数组
        for i = 1:current_batch_size
            sample_idx = start_idx + i - 1;
            all_orders(sample_idx) = batch_orders(i);
        end
        
        % 计算并显示当前进度
        progress = min(100, end_idx * 100 / num_samples);
        bar_filled = floor(progress * progress_bar_length / 100);
        progress_bar = ['[', repmat('#', 1, bar_filled), repmat(' ', 1, progress_bar_length - bar_filled), ']'];
        fprintf('\r进度 (σ_Δ = %.1f): %s %.1f%% (%d/%d)', sigma_delta, progress_bar, progress, end_idx, num_samples);
    end
    
    % 确保最后显示100%进度
    progress_bar = ['[', repmat('#', 1, progress_bar_length), ']'];
    fprintf('\r进度 (σ_Δ = %.1f): %s 100.0%% (%d/%d)\n', sigma_delta, progress_bar, num_samples, num_samples);
    
    % 统计twisted和nontwisted状态的比例
    twisted_count = sum(all_orders >= order_threshold);
    nontwisted_count = num_samples - twisted_count;
    
    twisted_ratio = twisted_count / num_samples;
    nontwisted_ratio = nontwisted_count / num_samples;
    
    % 存储结果（小数形式）
    results(delta_idx).sigma_delta = sigma_delta;
    results(delta_idx).twisted = twisted_ratio;
    results(delta_idx).nontwisted = nontwisted_ratio;
    
    % 计算并显示此sigma_delta的运行时间
    sigma_delta_time = toc(sigma_delta_start_time);
    fprintf('σ_Δ = %.1f 计算完成，用时 %.2f 秒 (%.2f 分钟)，twisted比例: %.4f\n', sigma_delta, sigma_delta_time, sigma_delta_time/60, twisted_ratio);
end

% 计算并显示总运行时间
total_time = toc(total_start_time);
fprintf('\n所有模拟完成，总运行时间: %.2f秒 (%.2f分钟)\n', total_time, total_time/60);

%% 准备绘图数据
sigma_delta_range = [results.sigma_delta];
twisted_data = [results.twisted];
nontwisted_data = [results.nontwisted];

%% 设置图表通用样式
set(0, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontName', 'Arial');
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontSize', 14);
set(0, 'DefaultLineLineWidth', 2.0);

%% 绘制twisted和nontwisted占比随sigma_delta变化的函数图
figure('Position', [100, 100, 1000, 700], 'Color', 'white');

% 使用两种颜色：一种是新颜色用于twisted，灰色用于nontwisted
colors = [
    [255, 140, 0]/255;     % 橙色 (twisted) - 与之前列出的所有颜色不同
    [127, 127, 127]/255    % 灰色 (nontwisted)
];

% 绘制数据点
hold on;
h1 = scatter(sigma_delta_range, twisted_data, 80, colors(1,:), 'filled', 'DisplayName', 'twisted');
h2 = scatter(sigma_delta_range, nontwisted_data, 80, colors(2,:), 'filled', 'DisplayName', 'nontwisted');

% 拟合曲线部分 - 使用平滑插值
if length(sigma_delta_range) >= 3
    finer_sigma_delta = linspace(min(sigma_delta_range)-0.1, max(sigma_delta_range)+0.1, 5000); % 增加采样点数量，使曲线更平滑
    
    % 使用平滑插值
    twisted_fit = smooth_interpolation(sigma_delta_range, twisted_data, finer_sigma_delta);
    nontwisted_fit = smooth_interpolation(sigma_delta_range, nontwisted_data, finer_sigma_delta);
    
    % 使用较宽的线条增强视觉效果
    plot(finer_sigma_delta, twisted_fit, 'Color', colors(1,:), 'LineWidth', 3.5, 'HandleVisibility', 'off');
    plot(finer_sigma_delta, nontwisted_fit, 'Color', colors(2,:), 'LineWidth', 3.5, 'HandleVisibility', 'off');
else
    % 直接连线部分
    plot(sigma_delta_range, twisted_data, 'Color', colors(1,:), 'LineWidth', 3.5, 'HandleVisibility', 'off');
    plot(sigma_delta_range, nontwisted_data, 'Color', colors(2,:), 'LineWidth', 3.5, 'HandleVisibility', 'off');
end
hold off;

% 添加网格和图例
grid on;
legend('Location', 'eastoutside', 'FontSize', 12);
xlabel('Three-Body Coupling Strength (σ_{\Delta})', 'FontSize', 16);
ylabel('Ratio', 'FontSize', 16); % 修改为比例而非百分比
title('Distribution of Twisted vs. Non-twisted States', 'FontSize', 20, 'FontWeight', 'bold');
subtitle(sprintf('No noise, t = %.0f sec, P_{order} threshold = %.2f', T_sim, order_threshold), 'FontSize', 14);

% 设置坐标轴范围和样式
xlim([min(sigma_delta_range) - 0.2, max(sigma_delta_range) + 0.2]);
ylim([-0.05, 1.05]); % 设置y轴范围从0到1

% 美化坐标轴
set(gca, 'Box', 'on', 'LineWidth', 1.5);
% 设置更精细的网格间隔
set(gca, 'XMinorTick', 'on', 'YMinorTick', 'on', 'XGrid', 'on', 'YGrid', 'on');
set(gca, 'GridLineStyle', '-', 'MinorGridLineStyle', ':');
set(gca, 'GridAlpha', 0.15, 'MinorGridAlpha', 0.1);

% 设置Y轴刻度
set(gca, 'YTick', 0:0.2:1);

% 优化图形渲染
set(gcf, 'Renderer', 'painters');  % 使用矢量渲染器以获得更清晰的线条

% 保存结果数据
save('twisted_vs_nontwisted_sigma_delta.mat', 'results', 'sigma_delta_range', 'twisted_data', 'nontwisted_data');

%% 输出表格格式的结果
fprintf('\n结果摘要（小数比例）:\n');
fprintf('σ_Δ\ttwisted\tnontwisted\t总比例\n');
fprintf('----------------------------------------\n');
for i = 1:length(results)
    fprintf('%.1f\t%.4f\t%.4f\t\t%.4f\n', ...
        results(i).sigma_delta, results(i).twisted, results(i).nontwisted, ...
        results(i).twisted + results(i).nontwisted);
end

%% 辅助函数：平滑插值
function y_interp = smooth_interpolation(x, y, x_interp)
    % 确保x和y是列向量
    x = x(:);
    y = y(:);
    
    % 使用MATLAB的内置插值和平滑功能
    try
        % 首先使用样条插值
        y_raw = interp1(x, y, x_interp, 'spline');
        
        % 然后对插值结果进行平滑处理
        window_size = max(5, round(length(y_raw) * 0.05)); % 平滑窗口大小，5或数据的5%
        y_interp = smooth(y_raw, window_size, 'loess'); % 使用局部加权回归平滑
    catch
        % 如果上述方法失败，回退到简单的三次样条插值
        warning('局部加权回归平滑失败，使用简单插值');
        y_interp = interp1(x, y, x_interp, 'pchip'); % 使用pchip插值（保形插值）
    end
    
    % 确保所有值都在[0,1]范围内
    y_interp = max(0, min(1, y_interp));
    
    % 确保输出尺寸与x_interp匹配
    y_interp = reshape(y_interp, size(x_interp));
end

%% 辅助函数：计算相互作用
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

%% 辅助函数：计算有序度 - 局部有序振荡器比例
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